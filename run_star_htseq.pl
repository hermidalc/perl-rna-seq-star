#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use sigtrap qw(handler sig_handler normal-signals error-signals ALRM);
use Cwd qw(cwd);
use Digest::MD5 qw(md5_hex);
use File::Basename qw(fileparse);
use File::Copy qw(move);
use File::Fetch;
use File::Path qw(make_path remove_tree);
use File::Spec;
use Getopt::Long qw(:config auto_help auto_version);
use IPC::Cmd qw(can_run);
use List::Util qw(max min notall uniq);
use Parallel::ForkManager;
use Pod::Usage qw(pod2usage);
use POSIX qw(SIGINT);
use Storable qw(lock_nstore lock_retrieve);
use Term::ANSIColor;
use Text::CSV qw(csv);
use Unix::Processors;
use Data::Dumper;

sub sig_handler {
    die "\nSTAR-HTSeq pipeline exited [", scalar localtime, "]\n";
}

our $VERSION = '0.1';

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Deepcopy = 1;

# unbuffer error and output streams (make sure STDOUT
# is last so that it remains the default filehandle)
select(STDERR); $| = 1;
select(STDOUT); $| = 1;

# config
my $sra_ftp_run_url_prefix =
    'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/';
my $fastq_dump_suffix_fmt = '_pass_';
my @bam_rg2sra_fields = (
    { RG => 'ID', SRA => 'Run' },
    { RG => 'SM', SRA => 'Sample' },
    { RG => 'LB', SRA => 'LibraryName' },
    { RG => 'PL', SRA => 'Platform' },
    { RG => 'PM', SRA => 'Model' },
);
my $procs = Unix::Processors->new();

# options
my $sra_query_str = '';
my $srr_file = '';
my $out_dir = File::Spec->abs2rel();
my $tmp_dir = cwd();
my $num_threads = $procs->max_physical;
my $genome_dir = 'star_grch38p2_d1_vd1_gtfv22';
my $gtf_file = 'gencode.v22.annotation.gtf';
my $use_sra_ftp = 0;
my $genome_shm = 0;
my $refresh_meta = 0;
my $query_only = 0;
my $regen_all = 0;
my $keep_all = 0;
my $keep_bams = 0;
my $gen_tx_bam = 0;
my $no_htseq = 0;
my $htseq_mode = 'intersection-nonempty';
my $htseq_stranded = 'no';
my $dry_run = 0;
my $verbose = 0;
my $debug = 0;
GetOptions(
    'sra-query:s' => \$sra_query_str,
    'srr-file:s' => \$srr_file,
    'out-dir:s' => \$out_dir,
    'tmp-dir:s' => \$tmp_dir,
    'num-threads:i' => \$num_threads,
    'genome-dir:s' => \$genome_dir,
    'gtf-file:s' => \$gtf_file,
    'use-sra-ftp' => \$use_sra_ftp,
    'genome-shm' => \$genome_shm,
    'refresh-meta' => \$refresh_meta,
    'query-only' => \$query_only,
    'regen-all' => \$regen_all,
    'keep-all' => \$keep_all,
    'keep-bams' => \$keep_bams,
    'gen-tx-bam' => \$gen_tx_bam,
    'no-htseq' => \$no_htseq,
    'htseq-mode' => \$htseq_mode,
    'htseq-stranded' => \$htseq_stranded,
    'dry-run' => \$dry_run,
    'verbose' => \$verbose,
    'debug' => \$debug,
) || pod2usage(-verbose => 0);
if (!$sra_query_str) {
    if (!$srr_file) {
        pod2usage(-message => 'Required --sra-query or --srr-file');
    }
    elsif (!-f $srr_file) {
        pod2usage(-message => 'Invalid --srr-file');
    }
}
if ($num_threads =~ /^\d+$/) {
    $num_threads = min(max(0, $num_threads), $procs->max_physical);
}
else {
    pod2usage(-message => '--num-threads must be an integer');
}
if (!$dry_run) {
    make_path($out_dir) unless -d $out_dir;
    make_path($tmp_dir) unless -d $tmp_dir;
}
$keep_bams = 0 if $keep_all;
print "#", '-' x 120, "#\n",
      "# STAR-HTSeq pipeline [" . scalar localtime() . "]\n\n";
my $srr_metadata;
if (!$sra_query_str) {
    print "Reading $srr_file: ";
    open(my $fh, '<', $srr_file);
    chomp(my @srr_ids = <$fh>);
    close($fh);
    @srr_ids = grep { /\S/ && /^SRR/ } @srr_ids;
    @srr_ids = sort { $a cmp $b } uniq @srr_ids;
    print scalar(@srr_ids), " unique SRRs\n";
    $sra_query_str = join(' OR ', map { "${_}[Accession]" } @srr_ids);
}
my $srr_metadata_pls_file =
    "$tmp_dir/".md5_hex($sra_query_str).".srr_metadata.pls";
if (!-f $srr_metadata_pls_file or $refresh_meta) {
    $sra_query_str =~ s/'/\'/g;
    my $sra_cmd_str = join(' | ',
        "esearch -db sra -query '$sra_query_str'",
        "efetch -format runinfo",
    );
    print 'Searching SRA:';
    print +($verbose or $debug) ? "\n$sra_cmd_str\n" : ' ';
    my $sra_meta_csv_str = `$sra_cmd_str`;
    die +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
        ": efetch failed (exit code ", $? >> 8, ")\n" if $?;
    $srr_metadata = csv(in => \$sra_meta_csv_str, headers => 'auto');
    my %seen_srrs;
    $srr_metadata = [ grep {
        $_->{'Run'} =~ /\S/ &&
        $_->{'Run'} =~ /^SRR/ &&
        !$seen_srrs{$_->{'Run'}}++
    } @{$srr_metadata} ];
    $srr_metadata = [
        sort { $a->{'Run'} cmp $b->{'Run'} } @{$srr_metadata}
    ];
    print scalar(@{$srr_metadata}), " runs found\n";
    print "Caching query metadata\n";
    if (!$dry_run) {
        lock_nstore($srr_metadata, $srr_metadata_pls_file);
    }
}
else {
    print "Loading cached query metadata\n";
    $srr_metadata = lock_retrieve($srr_metadata_pls_file);
    print scalar(@{$srr_metadata}), " runs found\n";
}
if ($debug) {
    print "SRR IDs: ", Dumper([ map { $_->{'Run'} } @{$srr_metadata} ]),
          "SRR metadata: ", Dumper($srr_metadata);
}
print "\n";
exit if $query_only;
my @srrs_completed = ();
my @htseq_cmd_data = ();
for my $srr_meta (@{$srr_metadata}) {
    my $srr_id = $srr_meta->{'Run'};
    print "[$srr_id]\n";
    my $srr_dir = File::Spec->abs2rel("$tmp_dir/$srr_id");
    my @fastq_files = map {
        "$srr_dir/${srr_id}${fastq_dump_suffix_fmt}${_}.fastq.gz"
    } 1..2;
    my $star_bam_file_name = 'Aligned.out.bam';
    my $star_bam_file = "$srr_dir/$star_bam_file_name";
    my $star_tx_bam_file_name = 'Aligned.toTranscriptome.out.bam';
    my $star_tx_bam_file = "$srr_dir/$star_tx_bam_file_name";
    my $star_counts_file_name = 'ReadsPerGene.out.tab';
    my $star_counts_file = "$srr_dir/$star_counts_file_name";
    my $star_completed_file = "$srr_dir/star_completed";
    my $final_star_bam_file_name = "${srr_id}.Aligned.bam";
    my $final_star_bam_file = "$out_dir/$final_star_bam_file_name";
    my $final_star_tx_bam_file_name = "${srr_id}.Aligned.toTranscriptome.bam";
    my $final_star_tx_bam_file = "$out_dir/$final_star_tx_bam_file_name";
    my $final_star_counts_file_name = "${srr_id}.star.counts.txt";
    my $final_star_counts_file = "$out_dir/$final_star_counts_file_name";
    my $htseq_counts_file_name = "${srr_id}.htseq.counts.txt";
    my $htseq_counts_file =
        ($keep_all ? $srr_dir : $out_dir)."/$htseq_counts_file_name";
    my $ran_star = 0;
    if (
        $regen_all or (
            (
                $keep_bams and (
                    (!-f $final_star_bam_file and !-f $star_bam_file) or
                    (
                        $gen_tx_bam and
                        !-f $final_star_tx_bam_file and !-f $star_tx_bam_file
                    )
                )
            ) or
            (!-f $final_star_counts_file and !-f $star_counts_file) or
            (-d $srr_dir and !-f $star_completed_file)
        )
    ) {
        if (!$dry_run) {
            if (-d $srr_dir) {
                remove_tree($srr_dir, { safe => 1 }) if $regen_all;
            }
            else {
                make_path($srr_dir);
            }
        }
        if ($regen_all or !-f "$srr_dir/${srr_id}.sra") {
            print "Downloading ${srr_id}.sra";
            my $dl_cmd_str;
            if ($use_sra_ftp) {
                my $sra_ftp_run_url = join('',
                    $sra_ftp_run_url_prefix, substr($srr_id, 0, 6),
                    "/$srr_id/${srr_id}.sra",
                );
                if (my $wget = can_run('wget')) {
                    $dl_cmd_str = "$wget -P '$srr_dir' $sra_ftp_run_url";
                }
                elsif (my $curl = can_run('curl')) {
                    $dl_cmd_str =
                        "$curl -o '$srr_dir/${srr_id}.sra' -L $sra_ftp_run_url";
                }
                else {
                    print "Using File::Fetch...\n";
                    my $ff = File::Fetch->new(uri => $sra_ftp_run_url);
                    if (!$dry_run) {
                        if (!$ff->fetch(to => $srr_dir)) {
                            exit(SIGINT) if ($? & 127) == SIGINT;
                            warn +(
                                -t STDERR ? colored('ERROR', 'red') : 'ERROR'
                            ), ": fetch failed, ", $ff->error, "\n\n";
                            next;
                        }
                    }
                }
            }
            else {
                $dl_cmd_str = "prefetch $srr_id -o '$srr_dir/${srr_id}.sra'";
            }
            if ($dl_cmd_str) {
                print "\n$dl_cmd_str" if $verbose or $debug;
                print "\n" unless $dl_cmd_str =~ /^prefetch /;
                if (!$dry_run) {
                    if (system($dl_cmd_str)) {
                        exit(SIGINT) if ($? & 127) == SIGINT;
                        warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                            ": download failed (exit code ", $? >> 8, ")\n\n";
                        next;
                    }
                }
                else {
                    print "\n";
                }
            }
        }
        else {
            print "Using existing ${srr_id}.sra\n";
        }
        if ($regen_all or notall { -f } @fastq_files) {
            my $val_cmd_str = join(' ', (
                "vdb-validate",
                "-x", "-B yes", "-I yes", "-C yes",
                "'$srr_dir/${srr_id}.sra'",
            ));
            print "Validating ${srr_id}.sra\n";
            print "$val_cmd_str\n" if $verbose or $debug;
            if (!$dry_run) {
                if (system($val_cmd_str)) {
                    exit(SIGINT) if ($? & 127) == SIGINT;
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                        ": vdb-validate failed (exit code ", $? >> 8, ")\n\n";
                    next;
                }
            }
            my $pfq_cmd_str = join(' ', (
                "parallel-fastq-dump",
                "--threads $num_threads",
                "--sra-id '$srr_dir/${srr_id}.sra'",
                "--outdir '$srr_dir'",
                "--tmpdir '$srr_dir'",
                "--gzip",
                "--skip-technical",
                "--readids",
                "--read-filter pass",
                "--dumpbase",
                "--split-3",
                "--clip",
            ));
            print "Generating ${srr_id}.sra FASTQs\n";
            print "$pfq_cmd_str\n" if $verbose or $debug;
            if (!$dry_run) {
                if (system($pfq_cmd_str)) {
                    exit(SIGINT) if ($? & 127) == SIGINT;
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                        ": parallel-fastq-dump failed (exit code ", $? >> 8,
                        ")\n\n";
                    next;
                }
            }
        }
        else {
            print "Using existing ",
                  join(' ', map { (fileparse($_))[0] } @fastq_files), "\n";
        }
        my @bam_rg_fields = map {
            "\"$_->{'RG'}:$srr_meta->{$_->{'SRA'}}\""
        } grep {
            $srr_meta->{$_->{'SRA'}} =~ /\S/
        } @bam_rg2sra_fields;
        my @quant_modes = ('GeneCounts');
        push @quant_modes, 'TranscriptomeSAM' if $gen_tx_bam;
        my $star_cmd_str = join(' ', (
            "STAR",
            "--runThreadN $num_threads",
            "--readFilesIn", join(' ', map { "'$_'" } @fastq_files),
            "--alignIntronMax 1000000",
            "--alignIntronMin 20",
            "--alignMatesGapMax 1000000",
            "--alignSJDBoverhangMin 1",
            "--alignSJoverhangMin 8",
            "--alignSoftClipAtReferenceEnds Yes",
            "--chimJunctionOverhangMin 15",
            "--chimMainSegmentMultNmax 1",
            "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip",
            "--chimSegmentMin 15",
            "--genomeDir '$genome_dir'",
            "--genomeLoad", $genome_shm ? "LoadAndKeep" : "NoSharedMemory",
            "--limitSjdbInsertNsj 1200000",
            "--outFileNamePrefix '$srr_dir/'",
            "--outSAMattrRGline", join(' ', @bam_rg_fields),
            "--outFilterIntronMotifs None",
            "--outFilterMatchNminOverLread 0.33",
            "--outFilterMismatchNmax 999",
            "--outFilterMismatchNoverLmax 0.1",
            "--outFilterMultimapNmax 20",
            "--outFilterScoreMinOverLread 0.33",
            "--outFilterType BySJout",
            "--outSAMattributes NH HI AS nM NM ch",
            "--outSAMstrandField intronMotif",
            "--outSAMtype BAM Unsorted",
            "--outSAMunmapped Within",
            "--quantMode", join(' ', @quant_modes),
            "--readFilesCommand zcat",
            "--twopassMode", $genome_shm ? "None" : "Basic",
        ));
        print "Running STAR alignment\n";
        print "$star_cmd_str\n" if $verbose or $debug;
        if (!$dry_run) {
            if (system($star_cmd_str)) {
                exit(SIGINT) if ($? & 127) == SIGINT;
                warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                    ": STAR failed (exit code ", $? >> 8, ")\n\n";
                next;
            }
            open(my $fh, '>', $star_completed_file);
        }
        $ran_star++;
    }
    if ($keep_bams) {
        if (
            $regen_all or !-f $final_star_bam_file or (
                -f $star_bam_file and
                (stat($final_star_bam_file))[7] !=
                (stat($star_bam_file))[7]
            )
        ) {
            if (!$ran_star) {
                print "Using existing $star_bam_file_name\n";
            }
            print "Moving $star_bam_file_name",
                  ' -> ',
                  "$final_star_bam_file_name\n";
            if (!$dry_run) {
                move($star_bam_file, $final_star_bam_file);
            }
        }
        else {
            print "Completed $final_star_bam_file_name exists\n";
        }
        if ($gen_tx_bam) {
            if (
                $regen_all or !-f $final_star_tx_bam_file or (
                    -f $star_tx_bam_file and
                    (stat($final_star_tx_bam_file))[7] !=
                    (stat($star_tx_bam_file))[7]
                )
            ) {
                if (!$ran_star) {
                    print "Using existing $star_tx_bam_file_name\n";
                }
                print "Moving $star_tx_bam_file_name",
                      ' -> ',
                      "$final_star_tx_bam_file_name\n";
                if (!$dry_run) {
                    move($star_tx_bam_file, $final_star_tx_bam_file);
                }
            }
            else {
                print "Completed $final_star_tx_bam_file_name exists\n";
            }
        }
    }
    if (
        $regen_all or !-f $final_star_counts_file or (
            -f $star_counts_file and
            (stat($final_star_counts_file))[7] !=
            (stat($star_counts_file))[7]
        )
    ) {
        if (!$ran_star) {
            print "Using existing $star_counts_file_name\n";
        }
        print "Moving $star_counts_file_name",
              ' -> ',
              "$final_star_counts_file_name\n";
        if (!$dry_run) {
            move($star_counts_file, $final_star_counts_file);
        }
    }
    else {
        print "Completed $final_star_counts_file_name exists\n";
    }
    if (!$no_htseq) {
        if (
            $regen_all or !-f $htseq_counts_file or
            (stat($htseq_counts_file))[7] == 0
        ) {
            if (!-f $star_bam_file) {
                $star_bam_file_name = $final_star_bam_file_name;
                $star_bam_file = $final_star_bam_file;
            }
            my $htseq_cmd_str = join(' ', (
                "htseq-count",
                "--format bam",
                "--order name",
                "--stranded $htseq_stranded",
                "--minaqual 10",
                "--type exon",
                "--idattr gene_id",
                "--mode $htseq_mode",
                $debug ? "" : "--quiet",
                "'$star_bam_file'",
                "'$gtf_file'",
            ));
            $htseq_cmd_str =~ s/\s+/ /g;
            if (!$keep_bams) {
                if (!$ran_star) {
                    print "Using existing $star_bam_file_name\n";
                }
                print "Running HTSeq quantification\n";
                print "$htseq_cmd_str\n" if $verbose or $debug;
                if (!$dry_run) {
                    open(my $fh, '>', $htseq_counts_file);
                    my $htseq_counts = `$htseq_cmd_str`;
                    if ($?) {
                        exit(SIGINT) if ($? & 127) == SIGINT;
                        warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                            ": htseq-count failed (exit code ", $? >> 8,
                            ")\n\n";
                        next;
                    }
                    else {
                        print $fh $htseq_counts;
                    }
                    close($fh);
                }
            }
            else {
                push @htseq_cmd_data, {
                    'srr_id' => $srr_id,
                    'cmd_str' => $htseq_cmd_str,
                    'counts_file' => $htseq_counts_file,
                }
            }
        }
        else {
            print "Completed $htseq_counts_file_name exists\n";
        }
    }
    if (!$dry_run and !$keep_all and -d $srr_dir) {
        remove_tree($srr_dir, { safe => 1 });
    }
    push @srrs_completed, $srr_id;
    print "\n";
}
if (@htseq_cmd_data) {
    print "Running HTSeq quantification\n";
    my $pm = Parallel::ForkManager->new($num_threads);
    for my $htseq_cmd (@htseq_cmd_data) {
        $pm->start and next;
        print "[$htseq_cmd->{'srr_id'}] Starting htseq-count\n";
        print "[$htseq_cmd->{'srr_id'}] $htseq_cmd->{'cmd_str'}\n"
            if $verbose or $debug;
        if (!$dry_run) {
            open(my $fh, '>', $htseq_cmd->{'counts_file'});
            my $htseq_counts = `$htseq_cmd->{'cmd_str'}`;
            if ($?) {
                exit(SIGINT) if ($? & 127) == SIGINT;
                warn "[$htseq_cmd->{'srr_id'}] ", (
                        -t STDERR ? colored('ERROR', 'red') : 'ERROR'
                     ), ": htseq-count failed (exit code ", $? >> 8, ")\n";
                next;
            }
            else {
                print $fh $htseq_counts;
            }
            close($fh);
        }
        print "[$htseq_cmd->{'srr_id'}] Finished htseq-count\n";
        $pm->finish;
    }
    $pm->wait_all_children;
    print "\n";
}
print scalar(@srrs_completed), " / ", scalar(@{$srr_metadata}),
      " SRRs completed\n\n",
      "STAR-HTSeq pipeline complete [", scalar localtime, "]\n\n";

__END__

=head1 NAME

run_star_htseq.pl - Run STAR-HTSeq Pipeline

=head1 SYNOPSIS

 run_star_htseq.pl [options]

 Options:
    --sra-query <str>    SRA query string to obtain SRR run metadata
                         (required if no --srr-file)
    --srr-file <file>    SRR ID list file
                         (required if no --sra-query)
    --out-dir <dir>      Output directory
                         (default = current dir)
    --tmp-dir <dir>      Temporary working directory
                         (default = current dir)
    --num-threads <n>    Number of parallel threads
                         (default = num cpus)
    --genome-dir <dir>   STAR genome index directory
                         (default = star_grch38p2_d1_vd1_gtfv22)
    --gtf-file <file>    Genome annotation GTF file
                         (default = gencode.v22.annotation.gtf)
    --use-sra-ftp        Use SRA direct FTP download instead of prefetch
                         (default = false)
    --genome-shm         Use STAR genome index in shared memory
                         (default = false)
    --refresh-meta       Re-query SRA to update metadata cache
                         (default = false)
    --query-only         Query SRA and then exit program
                         (default = false)
    --regen-all          Regenerate all result files
                         (default = false)
    --keep-all           Keep all result files
                         (default = false)
    --keep-bams          Keep STAR BAMs
                         (default = false)
    --gen-tx-bam         Generate STAR transcriptome-aligned BAM
                         (default = false)
    --no-htseq           Skip HTSeq quantification
                         (default = false)
    --htseq-mode         HTSeq --mode option
                         (default = intersection-nonempty)
    --htseq-stranded     HTSeq --stranded option
                         (default = no)
    --dry-run            Show what would've been done
                         (default = false)
    --verbose            Be verbose
    --help               Display usage and exit
    --version            Display program version and exit

=cut
