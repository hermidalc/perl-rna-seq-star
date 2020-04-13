# STAR + HTSeq + featureCounts RNA-seq Pipeline Wrapper

STAR + HTSeq + featureCounts RNA-seq processing pipeline environment and
wrapper script, including SRA query, download, and caching functionality and
useful reuse/restart features.

## Installation

Install latest
<a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">
Miniconda
</a>
for Python 3.

Clone the Github repository:

```bash
git clone git@github.com:hermidalc/rna-seq-star.git
```

Create `rna-seq-star` conda environment, either with Intel MKL:

```bash
cd rna-seq-star
conda env create -f envs/rna-seq-star-mkl.yml
conda activate rna-seq-star
```

or without:

```bash
cd rna-seq-star
conda env create -f envs/rna-seq-star.yml
conda activate rna-seq-star
```

Configure SRA toolkit:

The SRA toolkit prefetches files temporarily to the default cache directory
which might not be where you want.  If so, either set the `CACHE`
`public user-repository` location or under `TOOLS` prefetch
downloads to the `current directory`:

```bash
vdb-config -i
```

Install Aspera Connect (optional):

Downloading large .sra vdb files from SRA can sometimes result in repeated
download failures when using `prefetch` with the default download method. If
you install
<a href="https://www.ibm.com/aspera/connect/" target="_blank">
Aspera Connect
</a>
then `prefetch` will automatically use it and I've found it results in more
stable download performance.  To install Aspera Connect simply:

```bash
tar -xzf ibm-aspera-connect-version.tar.gz
./ibm-aspera-connect-version.sh
```

No more configuration is necessary as `prefetch` will automatically find and
use the `ascp` command.

## Usage

```
$ ./run_star.pl --help
Usage:
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
        --fcounts-stranded           featureCounts -s option
                                     (default = 0)
        --htseq                      Run HTSeq read quantification
                                     (default = true, false use --no-htseq)
        --htseq-par                  Run HTSeq jobs in parallel batches
                                     (default = true, false use --no-htseq-par)
        --htseq-num-par <n>          Number of HTSeq jobs in a batch
                                     (default = -1, num cpus)
        --htseq-mode                 HTSeq --mode option
                                     (default = intersection-nonempty)
        --htseq-stranded             HTSeq --stranded option
                                     (default = no)
        --min-aqual                  Minimum alignment quality score
                                     (default = 10)
        --dry-run                    Show what would've been done
                                     (default = false)
        --verbose                    Be verbose
        --help                       Display usage and exit
        --version                    Display program version and exit
```

The `--sra-query` can be any NCBI Entrez query string, e.g.

```bash
--sra-query 'SRP155030[All Fields] AND (1900[MDAT]:2900[MDAT] NOT "strategy exome"[Filter])'
```

#### Wrapper script features

The `run_star.pl` wrapper script encapsulates the following steps:

1.  Automatically creates a STAR genome index if doesn't already exist using the default or specified genome FASTA and GTF annotation file paths.

2.  Query SRA (or load SRR list) and then fetch and cache SRR metadata (with Entrez Direct tools)

Then for each SRR:

3.  Download SRA .sra file (via direct FTP URL, fallback to SRA toolkit prefetch) or with `--use-ena-fastqs` download ENA SRA FASTQs (via direct FTP URL) and skip to step 5
4.  Validate SRA vdb file (with SRA toolkit)
5.  Dump FASTQs from vdb (with parallel-fastq-dump wrapper to fastq-dump)
6.  Run STAR alignment and gene expression quantification
7.  Run featureCounts gene expression quantification
8.  Run HTSeq gene expression quantification

The wrapper script also has some useful features, for example:

*   Automatically caches downloaded SRA metadata from your query and will reuse cached metadata for the same SRA query.  Metadata can be refreshed with `--refresh-meta`.
*   Automatically looks for and reuses any existing and completed intermediate files from each processing step in case e.g. the script got killed. You simply restart the script with same parameters and it will come right back to where it left off. This can be overriden with `--regen-all`.
*   If any processing step fails for a run it will output the error details and automatically move on to the next run.
*   It saves a ton of space by only keeping the required result files necessary to complete all processing steps for each run.  This can be overriden by specifying file types to keep using  the `--keep ` option.
*   STAR detected novel splice junctions can be automatically filtered after the 1st pass with the option `--star-filter-sj-pass1`.
*   I've found HTSeq to be the slowest step in the pipeline by far, also because it is only single-threaded while the other compute intensive steps can be parallelized.  So to greatly improve performance the wrapper script by default runs HTSeq jobs in parallel batches after `--num-threads` runs have completed.  This can be overriden with `--no-htseq-par`.  You can make the batch size smaller than `--num-threads` with `--htseq-par-n`.
*   HTSeq quantification can be skipped with `--no-htseq` and then only STAR quantification files will be generated.

## References

<a href="https://www.ncbi.nlm.nih.gov/home/tools/" target="_blank">Entrez Direct utilities and SRA toolkit</a><br/>
<a href="https://github.com/rvalieris/parallel-fastq-dump" target="_blank">parallel-fastq-dump</a><br/>
<a href="https://github.com/alexdobin/STAR" target="_blank">STAR</a><br/>
<a href="https://github.com/simon-anders/htseq" target="_blank">HTSeq</a><br/>
<a href="http://subread.sourceforge.net/" target="_blank">Subread</a>
