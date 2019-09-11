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
use File::Find;
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
use URI;
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
my @states = qw(
    TMP_DIR SRA SRA_FASTQ ENA_FASTQ STAR MV_BAM MV_TX_BAM MV_ALL HTSEQ
);
my $sra_ftp_run_url_prefix =
    'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR';
my $ena_ftp_fastq_url_prefix = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq';
my @bam_rg2sra_fields = (
    { RG => 'ID', SRA => 'Run' },
    { RG => 'SM', SRA => 'Sample' },
    { RG => 'LB', SRA => 'LibraryName' },
    { RG => 'PL', SRA => 'Platform' },
    { RG => 'PM', SRA => 'Model' },
);
my $procs = Unix::Processors->new();

# options
my $sra_query = '';
my $srr_file = '';
my $out_dir = File::Spec->abs2rel();
my $tmp_dir = cwd();
my $num_threads = -1;
my $genome_dir = 'star_grch38p2_d1_vd1_gtfv22';
my $gtf_file = 'gencode.v22.annotation.gtf';
my $star_opts = '';
my @keep = ();
my $refresh_meta = 0;
my $query_only = 0;
my $use_ena_fastqs = 0;
my $genome_shm = 0;
my $regen_all = 0;
my $gen_tx_bam = 0;
my $htseq = 1;
my $htseq_par = 1;
my $htseq_mode = 'intersection-nonempty';
my $htseq_stranded = 'no';
my $dry_run = 0;
my $verbose = 0;
my $debug = 0;
GetOptions(
    'sra-query:s' => \$sra_query,
    'srr-file:s' => \$srr_file,
    'out-dir:s' => \$out_dir,
    'tmp-dir:s' => \$tmp_dir,
    'num-threads:i' => \$num_threads,
    'genome-dir:s' => \$genome_dir,
    'gtf-file:s' => \$gtf_file,
    'star-opts:s' => \$star_opts,
    'keep:s' => \@keep,
    'refresh-meta' => \$refresh_meta,
    'query-only' => \$query_only,
    'use-ena-fastqs' => \$use_ena_fastqs,
    'genome-shm' => \$genome_shm,
    'regen-all' => \$regen_all,
    'gen-tx-bam' => \$gen_tx_bam,
    'htseq!' => \$htseq,
    'htseq-par!' => \$htseq_par,
    'htseq-mode:s' => \$htseq_mode,
    'htseq-stranded:s' => \$htseq_stranded,
    'dry-run' => \$dry_run,
    'verbose' => \$verbose,
    'debug' => \$debug,
) || pod2usage(-verbose => 0);
if (!$sra_query) {
    if (!$srr_file) {
        pod2usage(-message => 'Required --sra-query or --srr-file');
    }
    elsif (!-f $srr_file) {
        pod2usage(-message => 'Invalid --srr-file');
    }
}
if ($num_threads =~ /^-?\d+$/) {
    $num_threads = $num_threads == -1
        ? $procs->max_physical
        : $num_threads > 0
            ? min($num_threads, $procs->max_physical)
            : $num_threads < -1
                ? max(1, $procs->max_physical + $num_threads + 1)
                : 1;
}
else {
    pod2usage(-message => '--num-threads must be an integer');
}
if (!$dry_run) {
    make_path($out_dir) unless -d $out_dir;
    make_path($tmp_dir) unless -d $tmp_dir;
}
@keep = split(' ', join(' ', @keep));
my %keep = map { $_ => 1 } @keep;
%keep = ( all => 1 ) if $keep{all};
print "#", '-' x 120, "#\n",
      "# STAR-HTSeq pipeline [" . scalar localtime() . "]\n\n";
my $srr_meta;
if (!$sra_query) {
    print "Reading $srr_file: ";
    open(my $fh, '<', $srr_file);
    chomp(my @srr_ids = <$fh>);
    close($fh);
    @srr_ids = grep { /\S/ && /^SRR/ } @srr_ids;
    @srr_ids = sort { $a cmp $b } uniq @srr_ids;
    print scalar(@srr_ids), " unique SRRs\n";
    $sra_query = join(' OR ', map { "$_\[Accession\]" } @srr_ids);
}
my $srr_meta_pls_file =
    "$tmp_dir/".md5_hex($sra_query).".srr_meta.pls";
if (!-f $srr_meta_pls_file or $refresh_meta) {
    $sra_query =~ s/'/\'/g;
    my $sra_cmd_str = join(' | ',
        "esearch -db sra -query '$sra_query'",
        "efetch -format runinfo",
    );
    print 'Searching SRA:';
    print +($verbose or $debug) ? "\n$sra_cmd_str\n" : ' ';
    my $sra_meta_csv_str = `$sra_cmd_str`;
    die +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
        ": efetch failed (exit code ", $? >> 8, ")\n" if $?;
    $srr_meta = csv(in => \$sra_meta_csv_str, headers => 'auto');
    my %seen_srrs;
    $srr_meta = [ grep {
        $_->{'Run'} =~ /\S/ &&
        $_->{'Run'} =~ /^SRR/ &&
        !$seen_srrs{$_->{'Run'}}++
    } @{$srr_meta} ];
    $srr_meta = [
        sort { $a->{'Run'} cmp $b->{'Run'} } @{$srr_meta}
    ];
    print scalar(@{$srr_meta}), " runs found\n";
    print "Caching query metadata\n";
    if (!$dry_run) {
        lock_nstore($srr_meta, $srr_meta_pls_file);
    }
}
else {
    print "Loading cached query metadata\n";
    $srr_meta = lock_retrieve($srr_meta_pls_file);
    print scalar(@{$srr_meta}), " runs found\n";
}
if ($debug) {
    print "SRR IDs: ", Dumper([ map { $_->{'Run'} } @{$srr_meta} ]),
          "SRR metadata: ", Dumper($srr_meta);
}
print "\n";
exit if $query_only;
my @srrs_completed = ();
my @htseq_run_data = ();
SRR: for my $run_idx (0 .. $#{$srr_meta}) {
    my $srr_id = $srr_meta->[$run_idx]->{'Run'};
    print "[$srr_id]\n";
    my $tmp_srr_dir = File::Spec->abs2rel("$tmp_dir/$srr_id");
    my $out_srr_dir = File::Spec->abs2rel("$out_dir/$srr_id");
    my $state_file = "$tmp_srr_dir/.state";
    my $init_state = (!$regen_all and -f $state_file)
        ? read_state($state_file) : {};
    my %tmp_file_name = (
        'sra' => "$srr_id.sra",
        (map {
            +"fastq$_" => join('',
                $srr_id,
                ($use_ena_fastqs and !$init_state->{SRA_FASTQ})
                    ? '_' : '_pass_',
                "$_.fastq.gz"
            )
        } 1..2),
        'star_bam' => 'Aligned.out.bam',
        'star_tx_bam' => 'Aligned.toTranscriptome.out.bam',
        'star_counts' => 'ReadsPerGene.out.tab',
        'htseq_counts' => 'htseq.counts.txt',
    );
    my %tmp_file = map {
        $_ => "$tmp_srr_dir/$tmp_file_name{$_}"
    } keys %tmp_file_name;
    my (%out_file_name, %out_file);
    if ($keep{all}) {
        %out_file_name = %tmp_file_name;
        %out_file = map {
            $_ => "$out_srr_dir/$out_file_name{$_}"
        } keys %out_file_name;
    }
    elsif ($keep{bam}) {
        %out_file_name = (
            'star_bam' => "${srr_id}.Aligned.bam",
            'star_tx_bam' => "${srr_id}.Aligned.toTranscriptome.bam",
            'star_counts' => "${srr_id}.star.counts.txt",
            'htseq_counts' => "${srr_id}.htseq.counts.txt",
        );
        %out_file = map {
            $_ => "$out_dir/$out_file_name{$_}"
        } keys %out_file_name;
        if (!$regen_all and !-f $state_file) {
            if (-f $out_file{'star_bam'} and
                (!$gen_tx_bam or -f $out_file{'star_tx_bam'}) and
                -f $out_file{'star_counts'}
            ) {
                $init_state = {
                    map { $_ => 1 } grep {
                        ($use_ena_fastqs
                            ? $_ !~ /^SRA(_FASTQ)?$/
                            : $_ ne 'ENA_FASTQ')
                        && $_ ne 'HTSEQ'
                    } @states
                };
                if ($htseq and -f $out_file{'htseq_counts'}) {
                    $init_state->{HTSEQ}++;
                }
            }
        }
    }
    else {
        %out_file_name = (
            'star_counts' => "${srr_id}.star.counts.txt",
            'htseq_counts' => "${srr_id}.htseq.counts.txt",
        );
        %out_file = map {
            $_ => "$out_dir/$out_file_name{$_}"
        } keys %out_file_name;
        for my $bam_type (qw(star_bam star_tx_bam)) {
            $out_file_name{$bam_type} = $tmp_file_name{$bam_type};
            $out_file{$bam_type} = $tmp_file{$bam_type};
        }
        if (!$regen_all and !-f $state_file) {
            if ($htseq) {
                if (
                    -f $out_file{'htseq_counts'} and
                    -f $out_file{'star_counts'}
                ) {
                    $init_state = {
                        map { $_ => 1 } grep {
                            ($use_ena_fastqs
                                ? $_ !~ /^SRA(_FASTQ)?$/
                                : $_ ne 'ENA_FASTQ')
                        } @states
                    };
                }
            }
            elsif (-f $out_file{'star_counts'}) {
                $init_state = {
                    map { $_ => 1 } grep {
                        ($use_ena_fastqs
                            ? $_ !~ /^SRA(_FASTQ)?$/
                            : $_ ne 'ENA_FASTQ')
                        && $_ ne 'HTSEQ'
                    } @states
                };
            }
        }
    }
    my $state = { %{$init_state} };
    if (!$state->{TMP_DIR}) {
        if (-d $tmp_srr_dir) {
            print "Cleaning directory $tmp_srr_dir\n";
            remove_tree($tmp_srr_dir, { safe => 1 }) unless $dry_run;
        }
        else {
            print "Creating directory $tmp_srr_dir\n";
            make_path($tmp_srr_dir) unless $dry_run;
        }
        $state->{TMP_DIR}++;
        write_state($state_file, $state) unless $dry_run;
    }
    elsif (!$state->{STAR} or !$state->{MV_ALL}) {
        print "Using existing directory $tmp_srr_dir\n";
    }
    if ($use_ena_fastqs and !$state->{SRA_FASTQ}) {
        if (!$state->{ENA_FASTQ}) {
            my $fastqs_downloaded = 0;
            for my $n (1..2) {
                print "Downloading $tmp_file_name{\"fastq$n\"}";
                if (download_url(
                    join('/',
                        $ena_ftp_fastq_url_prefix,
                        substr($srr_id, 0, 6),
                        '00'.substr($srr_id, -1),
                        $srr_id,
                        $tmp_file_name{"fastq$n"},
                    ),
                    $tmp_srr_dir,
                )) {
                    if (++$fastqs_downloaded == 2) {
                        $state->{ENA_FASTQ}++;
                        write_state($state_file, $state) unless $dry_run;
                    }
                }
            }
        }
        elsif (!$state->{STAR}) {
            print "Using existing ",
                  join(' ', @tmp_file_name{'fastq1', 'fastq2'}), "\n";
        }
    }
    if (!$state->{ENA_FASTQ}) {
        if (!$state->{SRA}) {
            print "Downloading $tmp_file_name{'sra'}";
            if (!download_url(
                join('/',
                    $sra_ftp_run_url_prefix,
                    substr($srr_id, 0, 6),
                    $srr_id,
                    "$tmp_file_name{'sra'}",
                ),
                $tmp_srr_dir,
            )) {
                my $dl_cmd_str =
                    "prefetch $srr_id -o '$tmp_file{'sra'}'";
                print "\n$dl_cmd_str" if $verbose or $debug;
                if (!$dry_run) {
                    if (system($dl_cmd_str)) {
                        exit($?) if ($? & 127) == SIGINT;
                        warn +(
                            -t STDERR ? colored('ERROR', 'red') : 'ERROR'
                        ), ": download failed (exit code ", $? >> 8, ")\n\n";
                        next SRR;
                    }
                }
                else {
                    print "\n";
                }
            }
            my $val_cmd_str = join(' ',
                "vdb-validate",
                "-x", "-B yes", "-I yes", "-C yes",
                "'$tmp_file{'sra'}'",
            );
            print "Validating $tmp_file_name{'sra'}\n";
            print "$val_cmd_str\n" if $verbose or $debug;
            if (!$dry_run) {
                if (system($val_cmd_str)) {
                    exit($?) if ($? & 127) == SIGINT;
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                        ": vdb-validate failed (exit code ", $? >> 8,
                        ")\n\n";
                    next SRR;
                }
            }
            $state->{SRA}++;
            write_state($state_file, $state) unless $dry_run;
        }
        elsif (!$state->{SRA_FASTQ}) {
            print "Using existing $tmp_file_name{'sra'}\n";
        }
        if (!$state->{SRA_FASTQ}) {
            my $pfq_cmd_str = join(' ',
                "parallel-fastq-dump",
                "--threads $num_threads",
                "--sra-id '$tmp_file{'sra'}'",
                "--outdir '$tmp_srr_dir'",
                "--tmpdir '$tmp_srr_dir'",
                "--gzip",
                "--skip-technical",
                "--readids",
                "--read-filter pass",
                "--dumpbase",
                "--split-3",
                "--clip",
            );
            print "Generating $tmp_file_name{'sra'} FASTQs\n";
            print "$pfq_cmd_str\n" if $verbose or $debug;
            if (!$dry_run) {
                if (system($pfq_cmd_str)) {
                    exit($?) if ($? & 127) == SIGINT;
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                        ": parallel-fastq-dump failed (exit code ",
                        $? >> 8, ")\n\n";
                    next SRR;
                }
            }
            $state->{SRA_FASTQ}++;
            write_state($state_file, $state) unless $dry_run;
        }
        elsif (!$state->{STAR}) {
            print "Using existing ",
                  join(' ', @tmp_file_name{'fastq1', 'fastq2'}), "\n";
        }
    }
    if (!$state->{STAR}) {
        my @bam_rg_fields = map {
            "\"$_->{'RG'}:$srr_meta->[$run_idx]->{$_->{'SRA'}}\""
        } grep {
            $srr_meta->[$run_idx]->{$_->{'SRA'}} =~ /\S/
        } @bam_rg2sra_fields;
        my @quant_modes = ('GeneCounts');
        push @quant_modes, 'TranscriptomeSAM' if $gen_tx_bam;
        my @star_cmd = (
            "STAR",
            "--runThreadN $num_threads",
            "--readFilesIn", join(' ',
                map { "'$_'" } @tmp_file{'fastq1', 'fastq2'}
            ),
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
            "--outFileNamePrefix '$tmp_srr_dir/'",
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
        );
        push @star_cmd, $star_opts if $star_opts;
        my $star_cmd_str = join(' ', @star_cmd);
        print "Running STAR alignment\n";
        print "$star_cmd_str\n" if $verbose or $debug;
        if (!$dry_run) {
            if (system($star_cmd_str)) {
                exit($?) if ($? & 127) == SIGINT;
                warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                    ": STAR failed (exit code ", $? >> 8, ")\n\n";
                next SRR;
            }
        }
        $state->{STAR}++;
        write_state($state_file, $state) unless $dry_run;
    }
    if ($keep{all}) {
        if (!$state->{MV_ALL}) {
            if ($out_srr_dir ne $tmp_srr_dir) {
                if (move_data($tmp_srr_dir, $out_srr_dir)) {
                    $state->{MV_ALL}++;
                    write_state($state_file, $state) unless $dry_run;
                }
                else {
                    next SRR;
                }
            }
        }
        else {
            print "Completed $out_srr_dir directory exists\n";
        }
    }
    elsif ($keep{bam}) {
        if (!$state->{MV_BAM}) {
            if (move_data($tmp_file{'star_bam'}, $out_file{'star_bam'})) {
                $state->{MV_BAM}++;
                write_state($state_file, $state) unless $dry_run;
            }
            else {
                next SRR;
            }
        }
        else {
            print "Completed $out_file_name{'star_bam'} exists\n";
        }
        if ($gen_tx_bam) {
            if (!$state->{MV_TX_BAM}) {
                if (move_data(
                    $tmp_file{'star_tx_bam'}, $out_file{'star_tx_bam'}
                )) {
                    $state->{MV_TX_BAM}++;
                    write_state($state_file, $state) unless $dry_run;
                }
                else {
                    next SRR;
                }
            }
            else {
                print "Completed $out_file_name{'star_tx_bam'} exists\n";
            }
        }
    }
    if (!$keep{all} or $keep{bam}) {
        if (!$state->{MV_ALL}) {
            if (move_data(
                $tmp_file{'star_counts'}, $out_file{'star_counts'}
            )) {
                $state->{MV_ALL}++;
                write_state($state_file, $state) unless $dry_run;
            }
            else {
                next SRR;
            }
        }
        else {
            print "Completed $out_file_name{'star_counts'} exists\n";
        }
    }
    if ($htseq) {
        if (!$state->{HTSEQ}) {
            my $htseq_cmd_str = join(' ',
                "htseq-count",
                "--format bam",
                "--order name",
                "--stranded $htseq_stranded",
                "--minaqual 10",
                "--type exon",
                "--idattr gene_id",
                "--mode $htseq_mode",
                $debug ? "" : "--quiet",
                "'$out_file{'star_bam'}'",
                "'$gtf_file'",
            );
            $htseq_cmd_str =~ s/\s+/ /g;
            if ($htseq_par) {
                push @htseq_run_data, {
                    'srr_id' => $srr_id,
                    'cmd_str' => $htseq_cmd_str,
                    'out_file_name' => $out_file_name{'htseq_counts'},
                    'out_file' => $out_file{'htseq_counts'},
                    'state' => $state,
                    'state_file' => $state_file,
                };
                if (
                    scalar(@htseq_run_data) % $num_threads == 0 or
                    $run_idx == $#{$srr_meta}
                ) {
                    print "\nRunning HTSeq quantification\n";
                    my $pm = Parallel::ForkManager->new($num_threads);
                    $pm->run_on_finish(sub {
                        my ($pid, $exit_code, $srr_id) = @_;
                        push @srrs_completed, $srr_id unless $exit_code;
                    });
                    HTSEQ: for my $htseq_run (@htseq_run_data) {
                        $pm->start($htseq_run->{'srr_id'}) and next HTSEQ;
                        print "[$htseq_run->{'srr_id'}] ",
                              "Starting htseq-count\n";
                        print "[$htseq_run->{'srr_id'}] ",
                              "$htseq_run->{'cmd_str'}\n"
                              if $verbose or $debug;
                        my $htseq_counts;
                        my $exit_code = 0;
                        if (!$dry_run) {
                            $htseq_counts = `$htseq_run->{'cmd_str'}`;
                            $exit_code = $?;
                        }
                        if ($exit_code) {
                            if (($exit_code & 127) == SIGINT) {
                                warn "htseq-count interrupted\n";
                            }
                            else {
                                warn +(-t STDERR
                                    ? colored('ERROR', 'red') :'ERROR'
                                ), ": htseq-count failed (exit code ",
                                    $exit_code >> 8, ")\n";
                            }
                        }
                        else {
                            print "[$htseq_run->{'srr_id'}] ",
                                  "Writing $htseq_run->{'out_file_name'}\n";
                            if (!$dry_run) {
                                open(my $fh, '>', $htseq_run->{'out_file'});
                                print $fh $htseq_counts;
                                close($fh);
                            }
                            $htseq_run->{'state'}->{HTSEQ}++;
                            if (!$dry_run and -f $htseq_run->{'state_file'}) {
                                write_state(
                                    $htseq_run->{'state_file'},
                                    $htseq_run->{'state'},
                                );
                            }
                        }
                        $pm->finish($exit_code);
                    }
                    $pm->wait_all_children;
                    @htseq_run_data = ();
                }
            }
            else {
                print "Using existing $out_file_name{'star_bam'}\n"
                    if $init_state->{STAR};
                print "Running HTSeq quantification\n";
                print "$htseq_cmd_str\n" if $verbose or $debug;
                my $htseq_counts;
                my $exit_code = 0;
                if (!$dry_run) {
                    $htseq_counts = `$htseq_cmd_str`;
                    $exit_code = $?;
                }
                if ($exit_code) {
                    exit($exit_code) if ($exit_code & 127) == SIGINT;
                    warn +(-t STDERR
                        ? colored('ERROR', 'red') :'ERROR'
                    ), ": htseq-count failed (exit code ",
                        $exit_code >> 8, ")\n";
                    next SRR;
                }
                else {
                    print "Writing $out_file_name{'htseq_counts'}\n";
                    if (!$dry_run) {
                        open(my $fh, '>', $out_file{'htseq_counts'});
                        print $fh $htseq_counts;
                        close($fh);
                    }
                    $state->{HTSEQ}++;
                    write_state($state_file, $state) unless $dry_run;
                    push @srrs_completed, $srr_id;
                }
            }
        }
        else {
            print "Completed $out_file_name{'htseq_counts'} exists\n";
        }
    }
    else {
        push @srrs_completed, $srr_id;
    }
    if (!$dry_run and !$keep{all} and -d $tmp_srr_dir) {
        if (
            $htseq and $htseq_par and !$keep{bam} and
            $run_idx != $#{$srr_meta}
        ) {
            finddepth({
                no_chdir => 1,
                wanted => sub {
                    if (-f) {
                        return if $_ eq $tmp_file{'star_bam'} or
                                  $_ eq $tmp_file{'star_tx_bam'} or
                                  $_ eq $state_file;
                        unlink $_;
                    }
                    elsif (-d and $_ ne $tmp_srr_dir) {
                        remove_tree($_, { safe => 1 });
                    }
                },
            }, $tmp_srr_dir);
        }
        else {
            remove_tree($tmp_srr_dir, { safe => 1 });
        }
    }
    print "\n";
}
print scalar(@srrs_completed), " / ", scalar(@{$srr_meta}),
      " SRRs completed\n\n",
      "STAR-HTSeq pipeline complete [", scalar localtime, "]\n\n";

sub read_state {
    my ($file) = @_;
    open(my $fh, '<', $file);
    my $state = { map { $_ => 1 } split(' ', <$fh>) };
    close($fh);
    return $state;
}

sub write_state {
    my ($file, $state) = @_;
    open(my $fh, '>', $file);
    print $fh join(' ', keys %{$state});
    close($fh);
}

sub download_url {
    my ($url, $dir) = @_;
    my $dl_cmd_str;
    my $name = (URI->new($url)->path_segments)[-1];
    if (my $wget = can_run('wget')) {
        $dl_cmd_str = "$wget -O '$dir/$name' '$url'";
    }
    elsif (my $curl = can_run('curl')) {
        $dl_cmd_str =
            "$curl -o '$dir/$name' -L '$url'";
    }
    print "\n";
    if ($dl_cmd_str) {
        print "$dl_cmd_str\n" if $verbose or $debug;
        if (!$dry_run) {
            if (system($dl_cmd_str)) {
                exit($?) if ($? & 127) == SIGINT;
                warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                    ": download failed (exit code ", $? >> 8, ")\n\n";
                return 0;
            }
        }
    }
    else {
        print "Using File::Fetch...\n";
        my $ff = File::Fetch->new(uri => $url);
        if (!$dry_run) {
            if (!$ff->fetch(to => $dir)) {
                exit($?) if ($? & 127) == SIGINT;
                warn +(
                    -t STDERR ? colored('ERROR', 'red') : 'ERROR'
                ), ": fetch failed, ", $ff->error, "\n\n";
                return 0;
            }
        }
    }
    return 1;
}

sub move_data {
    my ($src, $dest) = @_;
    print "Moving $src -> $dest\n";
    if (!$dry_run) {
        if (!move($src, $dest)) {
            exit($?) if ($? & 127) == SIGINT;
            warn +(
                -t STDERR ? colored('ERROR', 'red') : 'ERROR'
            ), ": move failed: $!\n\n";
            return 0;
        }
    }
    return 1;
}

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
                         (default = -1 which means all cpus)
    --genome-dir <dir>   STAR genome index directory
                         (default = star_grch38p2_d1_vd1_gtfv22)
    --gtf-file <file>    Genome annotation GTF file
                         (default = gencode.v22.annotation.gtf)
    --star-opts <str>    Additional STAR options (quoted string)
                         (default = none)
    --keep <str>         Additional file types to keep (quoted string)
                         (default = none, possible: all sra fastq bam)
    --refresh-meta       Re-query SRA to update metadata cache
                         (default = false)
    --query-only         Query SRA and cache metadata then exit
                         (default = false)
    --use-ena-fastqs     Download ENA SRA FASTQs (with SRA fallback)
                         (default = false)
    --genome-shm         Use STAR genome index in shared memory
                         (default = false)
    --regen-all          Regenerate all result files
                         (default = false)
    --gen-tx-bam         Generate STAR transcriptome-aligned BAM
                         (default = false)
    --htseq              Run HTSeq read quantification
                         (default = true, false use --no-htseq)
    --htseq-par          Run HTSeq in parallel batches
                         (default = true, false use --no-htseq-par)
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
