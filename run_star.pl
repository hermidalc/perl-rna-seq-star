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
use LWP::UserAgent;
use Parallel::ForkManager;
use Pod::Usage qw(pod2usage);
use POSIX qw(SIGINT);
use Storable qw(lock_nstore lock_retrieve);
use Term::ANSIColor;
use Text::CSV qw(csv);
use Unix::Processors;
use URI::Split qw(uri_split uri_join);
use Data::Dumper;

sub sig_handler {
    die "\nSTAR pipeline exited [", scalar localtime, "]\n";
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
    TMP_DIR SRA SRA_FASTQ ENA_FASTQ STAR_PASS1 STAR_PASS2
    MV_BAM MV_TX_BAM MV_COUNTS MV_ALL FCOUNTS HTSEQ
);
my $sra_run_ftp_url_prefix =
    'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra';
my $ena_fastq_ftp_search_url =
    'https://www.ebi.ac.uk/ena/data/warehouse/filereport'.
    '?result=read_run&fields=fastq_ftp&accession=';
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
my $run_file = '';
my @run_ids = ();
my $out_dir = File::Spec->abs2rel();
my $tmp_dir = cwd();
my $num_threads = -1;
my @genome_fasta_files = ('GRCh38.d1.vd1.fa');
my $genome_dir = 'star_index_grch38p2_d1_vd1_gtfv22';
my $gtf_file = 'gencode.v22.annotation.gtf';
my $star_genome_opts = '';
my $star_opts = '';
my $star_max_readlen = 100;
my $star_filter_sj_pass1 = 0;
my @keep = ();
my $refresh_meta = 0;
my $query_only = 0;
my $use_ena_fastqs = 0;
my $genome_shm = 0;
my $regen_all = 0;
my $gen_tx_bam = 0;
my $fcounts = 1;
my $htseq = 1;
my $htseq_par = 1;
my $htseq_num_par = -1;
my $htseq_mode = 'intersection-nonempty';
my $stranded = 'no';
my $min_aqual = 10;
my $dry_run = 0;
my $verbose = 0;
my $debug = 0;
GetOptions(
    'sra-query:s' => \$sra_query,
    'run-file:s' => \$run_file,
    'run-ids:s' => \@run_ids,
    'out-dir:s' => \$out_dir,
    'tmp-dir:s' => \$tmp_dir,
    'num-threads:i' => \$num_threads,
    'genome-fasta-file:s' => \@genome_fasta_files,
    'genome-dir:s' => \$genome_dir,
    'gtf-file:s' => \$gtf_file,
    'star-genome-opts:s' => \$star_genome_opts,
    'star-opts:s' => \$star_opts,
    'star-max-readlen:i' => \$star_max_readlen,
    'star-filter-sj-pass1' => \$star_filter_sj_pass1,
    'keep:s' => \@keep,
    'refresh-meta' => \$refresh_meta,
    'query-only' => \$query_only,
    'use-ena-fastqs' => \$use_ena_fastqs,
    'genome-shm' => \$genome_shm,
    'regen-all' => \$regen_all,
    'gen-tx-bam' => \$gen_tx_bam,
    'fcounts!' => $fcounts,
    'htseq!' => \$htseq,
    'htseq-par!' => \$htseq_par,
    'htseq-num-par:i' => \$htseq_num_par,
    'htseq-mode:s' => \$htseq_mode,
    'stranded:s' => \$stranded,
    'min-aqual:i' => \$min_aqual,
    'dry-run' => \$dry_run,
    'verbose' => \$verbose,
    'debug' => \$debug,
) || pod2usage(-verbose => 0);
@run_ids = split(' ', join(' ', @run_ids));
@keep = split(' ', join(' ', @keep));
my %keep = map { $_ => 1 } @keep;
%keep = ( all => 1 ) if $keep{all};
if (!$sra_query) {
    if (!$run_file) {
        if (!@run_ids) {
            pod2usage(
                -message => 'One or more of required: '.
                    '--sra-query, --run-file or --run-ids'
            );
        }
    }
    elsif (!-f $run_file) {
        pod2usage(-message => 'Invalid --run-file');
    }
}
for my $genome_fasta_file (@genome_fasta_files) {
    if (!-f $genome_fasta_file) {
        pod2usage(
            -message => "Invalid --genome-fasta-file $genome_fasta_file"
        );
    }
}
if (!-f $gtf_file) {
    pod2usage(-message => 'Invalid --gtf-file');
}
if ($star_max_readlen < 1) {
    pod2usage(-message => '--star-max-readlen must be a positive integer');
}
if ($stranded !~ /^(yes|no|reverse)$/) {
    pod2usage(-message => '--stranded must be one of yes|no|reverse');
}
$num_threads = $num_threads == -1
    ? $procs->max_physical
    : $num_threads > 0
        ? min($num_threads, $procs->max_physical)
        : $num_threads < -1
            ? max(1, $procs->max_physical + $num_threads + 1)
            : 1;
$htseq_num_par = $htseq_num_par == -1
    ? $num_threads
    : $htseq_num_par > 0
        ? $htseq_num_par
        : $htseq_num_par < -1
            ? max(1, $num_threads + $htseq_num_par + 1)
            : 1;
my $ua = LWP::UserAgent->new();
print "#", '-' x 120, "#\n",
      "# STAR pipeline [" . scalar localtime() . "]\n\n";
if ($run_file or @run_ids) {
    print "Command line: ", scalar(@run_ids), " runs\n" if @run_ids;
    if ($run_file) {
        print "Reading $run_file: ";
        open(my $fh, '<', $run_file);
        my @file_run_ids = <$fh>;
        close($fh);
        print scalar(@file_run_ids), " runs\n";
        push @run_ids, @file_run_ids;
    }
    @run_ids = uniq sort { $a cmp $b }
        grep { /\S/ && /^[A-Z]RR/ }
        map { s/\s+//gr } @run_ids;
    print scalar(@run_ids), " total unique runs\n";
    my @sra_query_parts = ($sra_query) if $sra_query;
    push @sra_query_parts, map { "$_\[Accession\]" } @run_ids;
    $sra_query = join(' OR ', @sra_query_parts);
}
my $run_meta;
my $run_meta_pls_file = "$tmp_dir/".md5_hex($sra_query).".run_meta.pls";
if (!-f $run_meta_pls_file or $refresh_meta) {
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
    $run_meta = csv(in => \$sra_meta_csv_str, headers => 'auto');
    my %seen_runs;
    $run_meta = [ grep {
        $_->{'Run'} =~ /\S/ &&
        $_->{'Run'} =~ /^[A-Z]RR/ &&
        !$seen_runs{$_->{'Run'}}++
    } @{$run_meta} ];
    $run_meta = [
        sort { $a->{'Run'} cmp $b->{'Run'} } @{$run_meta}
    ];
    print scalar(@{$run_meta}), " runs found\n";
    print "Caching query metadata\n";
    if (!$dry_run) {
        lock_nstore($run_meta, $run_meta_pls_file);
    }
}
else {
    print "Loading cached query metadata\n";
    $run_meta = lock_retrieve($run_meta_pls_file);
    print scalar(@{$run_meta}), " runs found\n";
}
if ($debug) {
    print "Run IDs: ", Dumper([ map { $_->{'Run'} } @{$run_meta} ]),
          "Run metadata: ", Dumper($run_meta);
}
print "\n";
exit if $query_only;
if (!-d $genome_dir) {
    print "Creating STAR genome index $genome_dir\n";
    make_path($genome_dir);
    my @star_cmd = (
        "STAR",
        "--runThreadN $num_threads",
        "--runMode genomeGenerate",
        "--genomeDir '$genome_dir'",
        "--genomeFastaFiles", join(' ',
            map { "'$_'" } @genome_fasta_files
        ),
    );
    push @star_cmd, $star_genome_opts if $star_genome_opts;
    my $star_cmd_str = join(' ', @star_cmd);
    print "$star_cmd_str\n" if $verbose or $debug;
    if (!$dry_run) {
        if (system($star_cmd_str)) {
            exit($?) if ($? & 127) == SIGINT;
            die +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                ": STAR failed (exit code ", $? >> 8, ")\n\n";
        }
    }
    print "\n";
}
if (!$dry_run) {
    make_path($out_dir) unless -d $out_dir;
    make_path($tmp_dir) unless -d $tmp_dir;
}
my @runs_completed = ();
my @htseq_run_data = ();
RUN: for my $run_idx (0..$#{$run_meta}) {
    my $run_id = $run_meta->[$run_idx]->{'Run'};
    print "[$run_id]\n";
    my $tmp_run_dir = File::Spec->abs2rel("$tmp_dir/$run_id");
    my $out_run_dir = File::Spec->abs2rel("$out_dir/$run_id");
    my $state_file = "$tmp_run_dir/.state";
    my $init_state = (!$regen_all and -f $state_file)
        ? read_state($state_file) : {};
    my %tmp_file_name = (
        'sra' => "$run_id.sra",
        (map {
            +"fastq$_" => join('',
                $run_id,
                ($use_ena_fastqs and !$init_state->{SRA_FASTQ})
                    ? '_' : '_pass_',
                "$_.fastq.gz"
            )
        } 1..2),
        'star_bam' => 'Aligned.out.bam',
        'star_tx_bam' => 'Aligned.toTranscriptome.out.bam',
        'star_counts' => 'ReadsPerGene.out.tab',
        'subread_counts' => 'subread.counts.tsv',
        'htseq_counts' => 'htseq.counts.tsv',
    );
    my %tmp_file = map {
        $_ => "$tmp_run_dir/$tmp_file_name{$_}"
    } keys %tmp_file_name;
    my (%out_file_name, %out_file);
    if ($keep{all}) {
        %out_file_name = %tmp_file_name;
        %out_file = map {
            $_ => "$out_run_dir/$out_file_name{$_}"
        } keys %out_file_name;
    }
    elsif ($keep{bam}) {
        %out_file_name = (
            'star_bam' => "${run_id}.Aligned.bam",
            'star_tx_bam' => "${run_id}.Aligned.toTranscriptome.bam",
            'star_counts' => "${run_id}.star.counts.tsv",
            'subread_counts' => "${run_id}.subread.counts.tsv",
            'htseq_counts' => "${run_id}.htseq.counts.tsv",
        );
        %out_file = map {
            $_ => "$out_dir/$out_file_name{$_}"
        } keys %out_file_name;
        if (!$regen_all and !-f $state_file) {
            if (
                -f $out_file{'star_bam'} and
                (!$gen_tx_bam or -f $out_file{'star_tx_bam'}) and
                -f $out_file{'star_counts'}
            ) {
                $init_state = {
                    map { $_ => 1 } grep {
                        $_ ne 'TMP_DIR' &&
                        ($use_ena_fastqs
                            ? $_ !~ /^SRA(_FASTQ)?$/
                            : $_ ne 'ENA_FASTQ') &&
                        $_ ne 'FCOUNTS' &&
                        $_ ne 'HTSEQ'
                    } @states
                };
                if ($fcounts and -f $out_file{'subread_counts'}) {
                    $init_state->{FCOUNTS}++;
                }
                if ($htseq and -f $out_file{'htseq_counts'}) {
                    $init_state->{HTSEQ}++;
                }
            }
        }
    }
    else {
        %out_file_name = (
            'star_counts' => "${run_id}.star.counts.tsv",
            'subread_counts' => "${run_id}.subread.counts.tsv",
            'htseq_counts' => "${run_id}.htseq.counts.tsv",
        );
        %out_file = map {
            $_ => "$out_dir/$out_file_name{$_}"
        } keys %out_file_name;
        for my $bam_type (qw(star_bam star_tx_bam)) {
            $out_file_name{$bam_type} = $tmp_file_name{$bam_type};
            $out_file{$bam_type} = $tmp_file{$bam_type};
        }
        if (!$regen_all and !-f $state_file) {
            if (-f $out_file{'star_counts'}) {
                $init_state = {
                    map { $_ => 1 } grep {
                        $_ ne 'TMP_DIR' &&
                        ($use_ena_fastqs
                            ? $_ !~ /^SRA(_FASTQ)?$/
                            : $_ ne 'ENA_FASTQ') &&
                        $_ ne 'FCOUNTS' &&
                        $_ ne 'HTSEQ'
                    } @states
                };
                if ($fcounts and -f $out_file{'subread_counts'}) {
                    $init_state->{FCOUNTS}++;
                }
                if ($htseq and -f $out_file{'htseq_counts'}) {
                    $init_state->{HTSEQ}++;
                }
            }
        }
    }
    my $state = { %{$init_state} };
    if (!$state->{STAR_PASS2}) {
        if (!$state->{TMP_DIR}) {
            if (-d $tmp_run_dir) {
                print "Cleaning directory $tmp_run_dir\n";
                remove_tree($tmp_run_dir, { safe => 1 }) unless $dry_run;
            }
            else {
                print "Creating directory $tmp_run_dir\n";
                make_path($tmp_run_dir) unless $dry_run;
            }
            $state->{TMP_DIR}++;
            write_state($state_file, $state) unless $dry_run;
        }
        else {
            print "Using existing directory $tmp_run_dir\n";
        }
        if ($use_ena_fastqs and !$state->{SRA_FASTQ}) {
            if (!$state->{ENA_FASTQ}) {
                print 'Querying ENA: ';
                my $response = $ua->get($ena_fastq_ftp_search_url . $run_id);
                if ($response->is_success) {
                    if ($response->decoded_content ne '') {
                        print "OK\n";
                        my @response_lines = split(
                            "\n", $response->decoded_content
                        );
                        my @fastq_ftp_urls = map { "ftp://$_" } split(
                            ';', $response_lines[$#response_lines], 2
                        );
                        my $fastqs_downloaded = 0;
                        for my $n (1..2) {
                            print "Downloading $tmp_file_name{\"fastq$n\"}\n";
                            if (download_url(
                                $fastq_ftp_urls[$n - 1],
                                $tmp_run_dir,
                                $tmp_file_name{"fastq$n"},
                            )) {
                                if (++$fastqs_downloaded == 2) {
                                    $state->{ENA_FASTQ}++;
                                    write_state($state_file, $state)
                                        unless $dry_run;
                                }
                            }
                            elsif (-f $tmp_file{"fastq$n"})  {
                                print "Removing $tmp_file_name{\"fastq$n\"}\n";
                                unlink $tmp_file{"fastq$n"};
                            }
                        }
                    }
                    else {
                        warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                             ": no FASTQs found for $run_id\n";
                    }
                }
                else {
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                         ': failed', $response->status_line, "\n";
                }
            }
            else {
                print "Using existing ",
                      join(' ', @tmp_file_name{'fastq1', 'fastq2'}), "\n";
            }
        }
        if (!$state->{ENA_FASTQ}) {
            if (!$state->{SRA}) {
                print "Downloading $tmp_file_name{'sra'}\n";
                if (!download_url(
                    join('/',
                        $sra_run_ftp_url_prefix,
                        substr($run_id, 0, 3),
                        substr($run_id, 0, 6),
                        $run_id,
                        $tmp_file_name{'sra'},
                    ),
                    $tmp_run_dir,
                    $tmp_file_name{'sra'},
                )) {
                    if (-f $tmp_file{'sra'}) {
                        print "Removing $tmp_file_name{'sra'}\n";
                        unlink $tmp_file{'sra'};
                    }
                    print 'Using SRA Toolkit prefetch';
                    my $dl_cmd_str =
                        "prefetch $run_id -o '$tmp_file{'sra'}'";
                    print "\n$dl_cmd_str" if $verbose or $debug;
                    if (!$dry_run) {
                        if (system($dl_cmd_str)) {
                            exit($?) if ($? & 127) == SIGINT;
                            warn +(
                                -t STDERR ? colored('ERROR', 'red') : 'ERROR'
                            ), ": download failed (exit code ", $? >> 8,
                            ")\n\n";
                            next RUN;
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
                        next RUN;
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
                    "--outdir '$tmp_run_dir'",
                    "--tmpdir '$tmp_run_dir'",
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
                        next RUN;
                    }
                }
                $state->{SRA_FASTQ}++;
                write_state($state_file, $state) unless $dry_run;
            }
            else {
                print "Using existing ",
                      join(' ', @tmp_file_name{'fastq1', 'fastq2'}), "\n";
            }
        }
        my $pass1_dir_name = '_STARpass1';
        my $pass1_dir = "$tmp_run_dir/$pass1_dir_name";
        my $pass1_sj_file_name = 'SJ.out.tab';
        my $pass1_sj_file = "$pass1_dir/$pass1_sj_file_name";
        my $pass1_sj_filter_file_name = 'SJ.filtered.out.tab';
        my $pass1_sj_filter_file = "$pass1_dir/$pass1_sj_filter_file_name";
        my @bam_rg_fields = map {
            "\"$_->{'RG'}:$run_meta->[$run_idx]->{$_->{'SRA'}}\""
        } grep {
            $run_meta->[$run_idx]->{$_->{'SRA'}} =~ /\S/
        } @bam_rg2sra_fields;
        my @quant_modes = ('GeneCounts');
        push @quant_modes, 'TranscriptomeSAM' if $gen_tx_bam;
        if ($star_filter_sj_pass1) {
            if (!$state->{STAR_PASS1}) {
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
                    "--genomeDir '$genome_dir'",
                    "--genomeLoad", $genome_shm
                        ? "LoadAndKeep" : "NoSharedMemory",
                    "--limitSjdbInsertNsj 1200000",
                    "--outFileNamePrefix '$pass1_dir/'",
                    "--outFilterIntronMotifs None",
                    "--outFilterMatchNminOverLread 0.33",
                    "--outFilterMismatchNmax 999",
                    "--outFilterMismatchNoverLmax 0.1",
                    "--outFilterMultimapNmax 20",
                    "--outFilterScoreMinOverLread 0.33",
                    "--outSAMtype None",
                    "--readFilesCommand zcat",
                    "--sjdbGTFfile '$gtf_file'",
                    "--sjdbOverhang", $star_max_readlen - 1,
                );
                my $star_cmd_str = join(' ', @star_cmd);
                print "Running STAR alignment pass 1\n";
                print "$star_cmd_str\n" if $verbose or $debug;
                if (!$dry_run) {
                    mkdir $pass1_dir, 0755;
                    if (system($star_cmd_str)) {
                        exit($?) if ($? & 127) == SIGINT;
                        warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                            ": STAR failed (exit code ", $? >> 8, ")\n\n";
                        next RUN;
                    }
                }
                print "Filtering $pass1_sj_file_name: ";
                my $num_novel_sj = 0;
                my $num_filter_novel_sj = 0;
                open(my $in_fh, '<', $pass1_sj_file);
                open(my $out_fh, '>', $pass1_sj_filter_file);
                while (<$in_fh>) {
                    chomp;
                    my @fields = split /\s+/;
                    # skip annotated
                    next unless $fields[5] == 0;
                    if (
                        # chromosomal and non-mitochondrial
                        $fields[0] =~ /chr([1-9][0-9]?|X|Y)/ and
                        # canonical
                        $fields[4] > 0 and
                        # supported by a unique mapper
                        $fields[6] > 0
                    ) {
                        print $out_fh "$_\n";
                        $num_filter_novel_sj++;
                    }
                    $num_novel_sj++;
                }
                close($in_fh);
                close($out_fh);
                print "$num_filter_novel_sj / $num_novel_sj",
                      " novel junctions kept\n";
                $state->{STAR_PASS1}++;
                write_state($state_file, $state) unless $dry_run;
            }
            else {
                print "Using existing $pass1_sj_filter_file_name\n";
            }
        }
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
            "--outFileNamePrefix '$tmp_run_dir/'",
            "--outFilterIntronMotifs None",
            "--outFilterMatchNminOverLread 0.33",
            "--outFilterMismatchNmax 999",
            "--outFilterMismatchNoverLmax 0.1",
            "--outFilterMultimapNmax 20",
            "--outFilterScoreMinOverLread 0.33",
            "--outFilterType BySJout",
            "--outSAMattrRGline", join(' ', @bam_rg_fields),
            "--outSAMattributes NH HI AS nM NM ch",
            "--outSAMstrandField intronMotif",
            "--outSAMtype BAM Unsorted",
            "--outSAMunmapped Within",
            "--quantMode", join(' ', @quant_modes),
            "--readFilesCommand zcat",
            "--sjdbGTFfile '$gtf_file'",
            "--sjdbOverhang", $star_max_readlen - 1,
        );
        push @star_cmd, $star_filter_sj_pass1
            ? "--sjdbFileChrStartEnd '$pass1_sj_filter_file'"
            : ("--twopassMode", $genome_shm ? "None" : "Basic");
        push @star_cmd, $star_opts if $star_opts;
        my $star_cmd_str = join(' ', @star_cmd);
        print "Running STAR alignment",
            $star_filter_sj_pass1 ? " pass 2\n" : "\n";
        print "$star_cmd_str\n" if $verbose or $debug;
        if (!$dry_run) {
            if (system($star_cmd_str)) {
                exit($?) if ($? & 127) == SIGINT;
                warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                    ": STAR failed (exit code ", $? >> 8, ")\n\n";
                next RUN;
            }
        }
        $state->{STAR_PASS2}++;
        write_state($state_file, $state) unless $dry_run;
    }
    if ($keep{all}) {
        if (!$state->{MV_ALL}) {
            if ($out_run_dir ne $tmp_run_dir) {
                if ($init_state->{TMP_DIR}) {
                    print "Using existing directory $tmp_run_dir\n";
                }
                if (move_data($tmp_run_dir, $out_run_dir)) {
                    $state->{MV_ALL}++;
                    write_state($state_file, $state) unless $dry_run;
                }
                else {
                    next RUN;
                }
            }
            elsif ($init_state->{STAR_PASS2}) {
                print "Completed $out_run_dir directory exists\n";
            }
        }
        else {
            print "Completed $out_run_dir directory exists\n";
        }
    }
    elsif ($keep{bam}) {
        if (!$state->{MV_BAM}) {
            if ($init_state->{TMP_DIR}) {
                print "Using existing directory $tmp_run_dir\n";
            }
            if (move_data($tmp_file{'star_bam'}, $out_file{'star_bam'})) {
                $state->{MV_BAM}++;
                write_state($state_file, $state) unless $dry_run;
            }
            else {
                next RUN;
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
                    next RUN;
                }
            }
            else {
                print "Completed $out_file_name{'star_tx_bam'} exists\n";
            }
        }
    }
    if (!$keep{all} or $keep{bam}) {
        if (!$state->{MV_COUNTS}) {
            if (!$keep{bam} and $init_state->{TMP_DIR}) {
                print "Using existing directory $tmp_run_dir\n";
            }
            if (move_data(
                $tmp_file{'star_counts'}, $out_file{'star_counts'}
            )) {
                $state->{MV_COUNTS}++;
                write_state($state_file, $state) unless $dry_run;
            }
            else {
                next RUN;
            }
        }
        else {
            print "Completed $out_file_name{'star_counts'} exists\n";
        }
    }
    if ($fcounts) {
        if (!$state->{FCOUNTS}) {
            if ($init_state->{STAR_PASS2}) {
                print "Using existing $out_file_name{'star_bam'}\n";
            }
            my $fcounts_cmd_str = join(' ',
                "featureCounts",
                "-T $num_threads",
                "-Q $min_aqual",
                "-p",
                "-s", $stranded eq 'no' ? 0 : $stranded eq 'yes' ? 1 : 2,
                "-a '$gtf_file'",
                "-o '$out_file{'subread_counts'}'",
                "'$out_file{'star_bam'}'",
            );
            $fcounts_cmd_str =~ s/\s+/ /g;
            print "Running featreCounts quantification\n";
            print "$fcounts_cmd_str\n" if $verbose or $debug;
            if (!$dry_run) {
                if (system($fcounts_cmd_str)) {
                    exit($?) if ($? & 127) == SIGINT;
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                        ": featureCounts failed (exit code ", $? >> 8, ")\n\n";
                    next RUN;
                }
            }
            $state->{FCOUNTS}++;
            write_state($state_file, $state) unless $dry_run;
        }
        else {
            print "Completed $out_file_name{'subread_counts'} exists\n";
        }
    }
    if ($htseq) {
        if (!$state->{HTSEQ}) {
            my $htseq_cmd_str = join(' ',
                "htseq-count",
                "--format bam",
                "--order name",
                "--stranded $stranded",
                "--minaqual $min_aqual",
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
                    'run_id' => $run_id,
                    'cmd_str' => $htseq_cmd_str,
                    'tmp_run_dir' => $tmp_run_dir,
                    'out_file_name' => $out_file_name{'htseq_counts'},
                    'out_file' => $out_file{'htseq_counts'},
                    'state' => $state,
                    'state_file' => $state_file,
                };
                if (!$keep{all} and !$keep{bam}) {
                    if ($state->{TMP_DIR}) {
                        print "Cleaning directory $tmp_run_dir except BAMs\n";
                    }
                    if (!$dry_run and -d $tmp_run_dir) {
                        finddepth({
                            no_chdir => 1,
                            wanted => sub {
                                if (-f) {
                                    return if $_ eq $tmp_file{'star_bam'} or
                                              $_ eq $tmp_file{'star_tx_bam'} or
                                              $_ eq $state_file;
                                    unlink $_;
                                }
                                elsif (-d and $_ ne $tmp_run_dir) {
                                    remove_tree($_, { safe => 1 });
                                }
                            },
                        }, $tmp_run_dir);
                    }
                }
                if (
                    scalar(@htseq_run_data) % $htseq_num_par == 0 or
                    $run_idx == $#{$run_meta}
                ) {
                    print "\nRunning HTSeq quantification\n";
                    my $pm = Parallel::ForkManager->new($htseq_num_par);
                    $pm->run_on_finish(sub {
                        my ($pid, $exit_code, $run_id) = @_;
                        push @runs_completed, $run_id unless $exit_code;
                    });
                    HTSEQ: for my $htseq_run (@htseq_run_data) {
                        $pm->start($htseq_run->{'run_id'}) and next HTSEQ;
                        print "[$htseq_run->{'run_id'}] ",
                              "Running htseq-count\n";
                        print "[$htseq_run->{'run_id'}] ",
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
                                warn "[$htseq_run->{'run_id'}] ",
                                     "htseq-count interrupted\n";
                            }
                            else {
                                warn "[$htseq_run->{'run_id'}] ",
                                     +(-t STDERR
                                         ? colored('ERROR', 'red') : 'ERROR'
                                     ), ": htseq-count failed (exit code ",
                                     $exit_code >> 8, ")\n";
                            }
                        }
                        else {
                            print "[$htseq_run->{'run_id'}] ",
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
                            if (!$keep{all}) {
                                if ($htseq_run->{'state'}->{TMP_DIR}) {
                                    print "[$htseq_run->{'run_id'}] ",
                                          "Removing directory ",
                                          "$htseq_run->{'tmp_run_dir'}\n";
                                }
                                if (
                                    !$dry_run and
                                    -d $htseq_run->{'tmp_run_dir'}
                                ) {
                                    remove_tree(
                                        $htseq_run->{'tmp_run_dir'},
                                        { safe => 1 }
                                    );
                                }
                            }
                        }
                        $pm->finish($exit_code);
                    }
                    $pm->wait_all_children;
                    @htseq_run_data = ();
                }
            }
            else {
                if ($init_state->{STAR_PASS2}) {
                    print "Using existing $out_file_name{'star_bam'}\n";
                }
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
                        ? colored('ERROR', 'red') : 'ERROR'
                    ), ": htseq-count failed (exit code ",
                        $exit_code >> 8, ")\n";
                    next RUN;
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
                }
            }
        }
        else {
            print "Completed $out_file_name{'htseq_counts'} exists\n";
        }
    }
    if (!$htseq or $state->{HTSEQ} or !$htseq_par) {
        push @runs_completed, $run_id;
        if (!$keep{all}) {
            if ($state->{TMP_DIR}) {
                print "Removing directory $tmp_run_dir\n";
            }
            if (!$dry_run and -d $tmp_run_dir) {
                remove_tree($tmp_run_dir, { safe => 1 });
            }
        }
    }
    print "\n";
}
print scalar(@runs_completed), " / ", scalar(@{$run_meta}),
      " runs completed\n\n",
      "STAR pipeline complete [", scalar localtime, "]\n\n";

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
    my ($url, $dir, $name) = @_;
    my @url_parts = (uri_split($url))[0..2];
    my @path_segments = split('/', $url_parts[-1]);
    my $base_path = join('/', @path_segments[0..$#path_segments-1]);
    my $base_url = uri_join(@url_parts[0..1], $base_path);
    print 'Checking FTP: ';
    if ($ua->head($base_url)->is_success) {
        print "OK\n";
        my $dl_cmd_str;
        if (my $wget = can_run('wget')) {
            $dl_cmd_str = "$wget -O '$dir/$name' '$url'";
        }
        elsif (my $curl = can_run('curl')) {
            $dl_cmd_str = "$curl -o '$dir/$name' -L '$url'";
        }
        if ($dl_cmd_str) {
            print "$dl_cmd_str\n" if $verbose or $debug;
            if (!$dry_run) {
                if (system($dl_cmd_str)) {
                    exit($?) if ($? & 127) == SIGINT;
                    warn +(-t STDERR ? colored('ERROR', 'red') : 'ERROR'),
                        ": download failed (exit code ", $? >> 8, ")\n";
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
                    ), ": fetch failed, ", $ff->error, "\n";
                    return 0;
                }
            }
        }
        return 1;
    }
    else {
        print "$url doesn't exist\n";
        return 0;
    }
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

run_star.pl - Run STAR Pipeline

=head1 SYNOPSIS

 run_star.pl [options]

 Options:
    --sra-query <str>            SRA query string to obtain run metadata
                                 (required if no --run-file or --run-ids)
    --run-file <file>            Run ID list file
                                 (required if no --run-query or --run-ids)
    --run-ids <str>              Run IDs (quoted string)
                                 (required if no --run-query or --run-file)
    --out-dir <dir>              Output directory
                                 (default = current dir)
    --tmp-dir <dir>              Temporary working directory
                                 (default = current dir)
    --num-threads <n>            Number of parallel threads
                                 (default = -1, num cpus)
    --genome-fasta-file <file>   STAR genome fasta file
                                 can specify option multiple times
                                 (default = GRCh38.d1.vd1.fa)
    --genome-dir <dir>           STAR genome index directory
                                 (default = star_genome_grch38p2_d1_vd1_gtfv22)
    --gtf-file <file>            Genome annotation GTF file
                                 (default = gencode.v22.annotation.gtf)
    --star-genome-opts <str>     Additional STAR genome options (quoted string)
                                 (default = none)
    --star-opts <str>            Additional STAR mapping options (quoted string)
                                 (default = none)
    --star-max-readlen <n>       STAR maximum read length
                                 (default = 100)
    --star-filter-sj-pass1       Filter STAR 1st pass novel splice junctions
                                 for chromosomal (non-mito), canonical, and
                                 unique mapper supported splice junctions
                                 (default = false)
    --keep <str>                 Additional file types to keep (quoted string)
                                 (default = none, possible: all sra fastq bam)
    --refresh-meta               Re-query SRA to update metadata cache
                                 (default = false)
    --query-only                 Query SRA and cache metadata then exit
                                 (default = false)
    --use-ena-fastqs             Download ENA SRA FASTQs (with SRA fallback)
                                 (default = false)
    --genome-shm                 Use STAR genome index in shared memory
                                 (default = false)
    --regen-all                  Regenerate all result files
                                 (default = false)
    --gen-tx-bam                 Generate STAR transcriptome-aligned BAM
                                 (default = false)
    --fcounts                    Run featureCounts read quantification
                                 (default = true, false use --no-fcounts)
    --htseq                      Run HTSeq read quantification
                                 (default = true, false use --no-htseq)
    --htseq-par                  Run HTSeq jobs in parallel batches
                                 (default = true, false use --no-htseq-par)
    --htseq-num-par <n>          Number of HTSeq jobs in a batch
                                 (default = -1, num cpus)
    --htseq-mode                 HTSeq --mode option
                                 (default = intersection-nonempty)
    --stranded                   Library prep strand specificity
                                 (default = no, yes|no|reverse)
    --min-aqual                  Minimum alignment quality score
                                 (default = 10)
    --dry-run                    Show what would've been done
                                 (default = false)
    --verbose                    Be verbose
    --help                       Display usage and exit
    --version                    Display program version and exit

=cut
