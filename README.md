# STAR-HTSeq RNA-seq Pipeline Wrapper

STAR-HTSeq RNA-seq processing pipeline environment and wrapper script,
including SRA query, download, and caching functionality and useful
reuse/restart features.

## Installation

Install latest
<a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">
Miniconda
</a>
or
<a href="https://www.anaconda.com/distribution/" target="_blank">
Anaconda
</a>
for Python 3.x. and make sure `conda` is updated.

Clone the Github repository:

```bash
git clone git@github.com:hermidalc/rna-seq-star-htseq.git
```

Create `rna-seq-star-htseq` conda environment:

```bash
cd rna-seq-star-htseq
conda env create -f environment.yml
conda activate rna-seq-star-htseq
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
<a href="https://downloads.asperasoft.com/connect2/" target="_blank">
Aspera Connect
</a>
then `prefetch` will automatically use it and I've found it results in more
stable download performance.  To install Aspera Connect simply:

```bash
tar -xzf ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.tar.gz
./ibm-aspera-connect-3.9.6.173386-linux-g2.12-64.sh
```

No more configuration is necessary as `prefetch` will automatically find and
use the `ascp` command.

## Usage

#### Creating a STAR genome index

For example, for GRCh38 using the current NCI GDC FASTA and GENCODE GTF
versions (which will require ~36GB RAM to map reads using this genome):

```bash
mkdir star_grch38p2_d1_vd1_gtfv22

STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir star_grch38p2_d1_vd1_gtfv22 \
--genomeFastaFiles GRCh38.d1.vd1.fa \
--sjdbGTFfile gencode.v22.annotation.gtf \
--sjdbOverhang 100
```

#### Wrapper script usage

```
$ run_star_htseq.pl -h
Usage:
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
        --use-sra-ftp        Use SRA FTP URL download instead of prefetch
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
```

The `--sra-query` can be any NCBI Entrez query string, e.g.

```bash
--sra-query 'SRP155030[All Fields] AND (1900[MDAT]:2900[MDAT] NOT "strategy exome"[Filter])'
```

#### Wrapper script features

The `run_star_htseq.pl` wrapper script encapsulates the following steps:

1.  Query SRA (or load SRR list) and then fetch and cache SRR metadata (with Entrez Direct tools)

Then for each SRR:

2.  Prefetch SRA vdb file (with SRA toolkit)
3.  Validate SRA vdb file (with SRA toolkit)
4.  Dump FASTQs from vdb (with parallel-fastq-dump wrapper to fastq-dump)
5.  Run STAR alignment and gene expression quantification
6.  Run HTSeq gene expression quantification

The wrapper script also has some useful features, for example:

*   It automatically caches downloaded SRA metadata for your query and will reuse cached metadata for the same SRA query.  Metadata can be refreshed with `--refresh-meta`.
*   It automatically looks for and reuses any existing completed intermediate files from each processing step in case e.g. the script got killed. You simply restart the script with same parameters and it will come right back to where it left off. This can be overriden with `--regen-all`.
*   If the processing fails for any SRR it will output that there was an error and automatically move on to the next.
*   It saves a ton of space by not keeping intermediate result files after finishing processing for each SRR.  This can be overriden with the options `--keep-all` or for just BAMs `--keep-bams`.
*   In my environment I've found HTSeq to be the slowest step in the pipeline by far, also because it is only single-threaded while the other compute intensive steps are parallelized.  So to greatly improve performance the wrapper script can run the HTSeq jobs in parallel at the end if you specify the option `--keep-bams`, though beware that this will take up a lot more space.

## References

<a href="https://www.ncbi.nlm.nih.gov/home/tools/" target="_blank">Entrez Direct utilities and SRA toolkit</a><br/>
<a href="https://github.com/rvalieris/parallel-fastq-dump" target="_blank">parallel-fastq-dump</a><br/>
<a href="https://github.com/alexdobin/STAR" target="_blank">STAR</a><br/>
<a href="https://github.com/simon-anders/htseq" target="_blank">HTSeq</a>
