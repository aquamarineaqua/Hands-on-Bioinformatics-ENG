**Chapter 1** RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification

**Section 1-1** STAR alignment

This section is adapted from *STAT 115 2021 Homework 1 Problem 1*

## 目录

  - [1. Installing STAR](#1-installing-star)
  - [2. Data Preparation](#2-data-preparation)
    - [Three Required Data Files](#three-required-data-files)
    - [Previewing Genome Information](#previewing-genome-information)
  - [3 Index Generation](#3-index-generation)
  - [4 RNA-Seq Data Alignment](#4-rna-seq-data-alignment)
    - [(1) Obtaining Yeast RNA-Seq Raw Data](#1-obtaining-yeast-rna-seq-raw-data)
    - [(2) Running STAR Alignment](#2-running-star-alignment)
  - [5 Reviewing RNA-Seq Alignment Results](#5-reviewing-rna-seq-alignment-results)
    - [(1) Checking Log Files for Statistics](#1-checking-log-files-for-statistics)
      - [Example Interpretation](#example-interpretation-1m_srr9336468_logfinalout)
    - [(2) Viewing BAM Files](#2-viewing-bam-files)
      - [① Using `Samtools` on Linux](#①-using-samtools-on-linux)
      - [② Using Python Library `pysam`](#②-using-python-library-pysam)

From this chapter onward, we will process and analyze High-Throughput Sequencing (HTS) data, specifically **RNA-Seq** data, covering alignment, quality control, and quantification. This section focuses on practical **alignment** using `STAR`, with the yeast genome (*Saccharomyces cerevisiae*) as the study model.

## 1. Installing STAR

We create a new Conda environment `star_env` and install `STAR`:

```bash
conda create -n star_env -c bioconda star
conda activate star_env
```

Verify the installation:

```bash
STAR --version
```

## 2. Data Preparation

### Three Required Data Files

We need three types of data files:

1. **Genome Reference File**:
   - Represents the complete DNA sequence, used as a reference for RNA-Seq alignment.
   - Typically sourced from public databases such as **Ensembl**, e.g., the human genome `GRCh38`.
   - File format: *FASTA*.
2. **Annotation File**:
   - Corresponds to the genome reference file.
   - Example: `GRCh38.gtf` for the human genome from **Ensembl**.
   - **GTF (Gene Transfer Format)** files describe gene and transcript coordinates and features.
3. **RNA-Seq Raw Data**:
   - Typically in *FASTQ* format, generated by sequencing platforms (e.g., Illumina).
   - Example: `sample_1.fastq.gz`, used for alignment against the reference genome with `STAR`.

---

**Downloading the Reference Genome**

**(1)** Create a directory to store the data:

```
mkdir -p ~/STAR_data/genome
cd ~/STAR_data/genome
```

**(2)** Download the **Saccharomyces cerevisiae** genome and annotation files from *Ensembl* and decompress them:

```bash
wget ftp://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-107/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.107.gtf.gz
gunzip *.gz
```

### Previewing Genome Information

We can use Python's `biopython` library to read the genome FASTA file and inspect its basic details. Additionally, tools such as `IGV` can visualize whole genomes.

**Installing** `biopython`:

```bash
pip install biopython
```

**Previewing Genome Information**:

```python
from Bio import SeqIO

# Define file path
file_path = "genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

# Read and print basic information
with open(file_path, "r") as file:
    records = list(SeqIO.parse(file, "fasta"))
    
    # Summary of sequences
    print(f"共有 {len(records)} 条序列")  # Total sequences
    for i, record in enumerate(records[:5]):  # View first 5 records
        print(f"\n序列 {i + 1}:")  # Sequence
        print(f"ID: {record.id}")
        print(f"描述: {record.description}")  # Description
        print(f"序列长度: {len(record.seq)}")  # Length
        print(f"前 100 个碱基: {record.seq[:100]}")  # First 100 bases
```

![Image](https://github.com/user-attachments/assets/ca385422-587a-4b5e-b4bb-8169723d1222)

```python
# Extract sequence lengths and IDs
seq_lengths = [len(record.seq) for record in records]
seq_ids = [record.id for record in records]

import matplotlib.pyplot as plt
# Generate bar plot
plt.figure(figsize=(10, 6))
plt.barh(seq_ids, seq_lengths, color='skyblue')
plt.xlabel("Sequence Length (bp)")
plt.ylabel("Sequence ID")
plt.title("Overview of Saccharomyces cerevisiae Genome Sequences")
plt.tight_layout()
plt.show()
```

![Image](https://github.com/user-attachments/assets/1f4276ef-78c1-4b3b-bc52-f4837d551790)

## 3 Index Generation

The main purposes of **indexing** are:

1. Accelerating the process of mapping reads to the reference genome, optimizing alignment performance.
2. Supporting spliced alignment across exon junctions.

We will use the reference genome and annotation file to generate the STAR genome index.

```bash
STAR --runMode genomeGenerate \
     --genomeDir ~/STAR_data/index \
     --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
     --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.107.gtf \
     --runThreadN 4 \
     --genomeChrBinNbits 10 \
     --genomeSAindexNbases 10
```

**Parameter Explanation:**

- `--genomeDir`: Directory to store the generated index files.
- `--genomeFastaFiles`: Path to the genome *FASTA* file.
- `--sjdbGTFfile`: Path to the annotation *GTF* file.
- `--runThreadN`: Number of parallel threads (adjust based on CPU cores).
- `--genomeChrBinNbits`: Default is 18; reducing it to 10 decreases memory usage.

**What Files Are Generated After Indexing?**

The generated index files include:

- **`Genome`** (binary-encoded genome sequence)
- **`SA`** (suffix array)
- **`SAindex`** (FM-Index)
- **`sjdbInfo.txt`** (splice junction database)

These index files will be loaded by **STAR** during the alignment process to speed up read mapping.

![Image](https://github.com/user-attachments/assets/16ef8827-80f6-44ee-885e-e3b8ef7e763d)

## 4 RNA-Seq Data Alignment

### (1) Obtaining Yeast RNA-Seq Raw Data

We will use a set of paired-end yeast RNA-Seq data for STAR testing.

> **Dataset Information:**
>
> - Data source: [https://bioinfogp.cnb.csic.es/files/samples/rnaseq/](https://bioinfogp.cnb.csic.es/files/samples/rnaseq/)
> - The dataset **"RNA-Seq_Sample_Files.zip"** (1.3 GB) contains:
>   - **Sequencing Data**: Multiple sets of *Illumina* paired-end sequencing *FASTQ* files, with three biological replicates per condition.
>   - **Genome Sequence**: The **Saccharomyces cerevisiae** genome sequence *FASTA* file (previously downloaded).
>   - **Gene Annotation**: A *GTF* file with gene coordinates (previously downloaded) and a tab-separated text file (`ann.txt`) containing gene symbols and descriptions.
>
> This dataset is suitable for running a simple RNA-Seq workflow, making it useful for educational and research purposes.
>
> ![Image](https://github.com/user-attachments/assets/03794abd-3094-4ff0-a6eb-84fb4634974c)

**1 Creating a Storage Directory**

```bash
mkdir -p ~/STAR_data/test_rnaseq
cd ~/STAR_data/test_rnaseq
```

**2 Preparing Data**

We will select three paired-end replicates as test samples, representing different experimental conditions (pH and CO2 levels):

1. `1M_SRR9336468_1.fastq.gz` (Read 1) & `1M_SRR9336468_2.fastq.gz` (Read 2)
2. `1M_SRR9336471_1.fastq.gz` (Read 1) & `1M_SRR9336471_2.fastq.gz` (Read 2)
3. `1M_SRR9336474_1.fastq.gz` (Read 1) & `1M_SRR9336474_2.fastq.gz` (Read 2)

**Note:**

- **Paired-End Sequencing**: Each fragment is sequenced from both ends, generating two paired reads stored in separate files, typically named **Read 1** and **Read 2** (or **File 1** and **File 2**). These represent **Forward Reads** and **Reverse Reads**, respectively.
- **SRR Identifiers**: Each sequencing run has a unique identifier in the **SRA (Sequence Read Archive)** database. For example, `SRR9336468` refers to a specific sequencing experiment.

> **SRA** (Sequence Read Archive) is a public repository maintained by **NCBI**, storing RNA-seq and DNA-seq data from sequencing platforms like *Illumina*, *PacBio*, and *Nanopore*.
>
> - Official website: [Home - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra)

We will copy these six `gz` files into the `test_rnaseq` directory and extract them.

```bash
gunzip *.gz
```



### (2) Running STAR Alignment

We will use `STAR` to align this RNA-Seq dataset.

**First Sample `SRR9336468`**

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn 1M_SRR9336468_1.fastq 1M_SRR9336468_2.fastq \
     --runThreadN 4 \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336468_ \
     --outSAMtype BAM SortedByCoordinate
```

**Parameter Explanation**

- `--genomeDir`: Directory of the previously generated yeast genome index.
- `--readFilesIn`: Input paired-end FASTQ files.
- `--runThreadN`: Number of parallel threads.
- `--outFileNamePrefix`: Prefix for output files.
- `--outSAMtype`: Output format as BAM file, sorted by coordinate.

**Commands for Second and Third Samples `SRR9336471,SRR9336474`**

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn 1M_SRR9336471_1.fastq 1M_SRR9336471_2.fastq \
     --runThreadN 4 \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336471_ \
     --outSAMtype BAM SortedByCoordinate
     
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn 1M_SRR9336474_1.fastq 1M_SRR9336474_2.fastq \
     --runThreadN 4 \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336474_ \
     --outSAMtype BAM SortedByCoordinate
```

**Output Files**

After successful execution, the following results will be generated in `~/STAR_results/`:

- **Alignment Results (BAM File)**: `1M_SRR9336468_Aligned.sortedByCoord.out.bam`
- **Alignment Log File**: `1M_SRR9336468_Log.final.out`



## 5 Reviewing RNA-Seq Alignment Results

### (1) Checking Log Files for Statistics

Run the following command to check the alignment statistics:

```bash
cat ~/STAR_results/1M_SRR9336468_Log.final.out
cat ~/STAR_results/1M_SRR9336471_Log.final.out
cat ~/STAR_results/1M_SRR9336474_Log.final.out
```

![Image](https://github.com/user-attachments/assets/48b6a520-85ab-4484-8f60-84d381f2a123)

![Image](https://github.com/user-attachments/assets/dddf9b21-5181-478a-86db-36738669e98d)

![Image](https://github.com/user-attachments/assets/6b7f9394-fa3e-40f7-8745-e979a7479f47)

#### Example Interpretation (`1M_SRR9336468_Log.final.out`):

- Number of input reads: `1,000,000`
- Average input read length: `300 bp`

**Unique Mapping (Unique Reads)**

> **Definition:** Unique reads align to a single, specific location in the reference genome without any other matches or with significantly lower scores at other locations.

- Uniquely mapped reads number: `881,977`
- Uniquely mapped reads percentage: `88.20%`

**Information of Splices**

Includes:

- Number of splices: Total – Total detected splicing events

- Number of splices: Annotated (sjdb) – Splicing events matching the annotations in the GTF file

- Number of splices: Non-canonical – Non-canonical splicing sites

> **Splices** are the junctions between **exons** and **introns** in a gene. The detected splices count represents exon-exon junctions identified during alignment.
> **Canonical splice sites** follow the common **GT-AG** rule (GT at the intron start, AG at the end).
> **Non-canonical splice sites** may have alternative motifs (e.g., AT-AC or GC-AG).

**Mismatch, Deletion and Insertion**

- Mismatch rate per base: `0.25%`
- Deletion rate per base: `0.02%`
- Insertion rate per base: `0.01%`

**Multi-Mapping Reads**

> **Definition:** Reads that align to multiple locations in the genome with similar alignment scores, typically originating from repetitive regions.

- Number of reads mapped to multiple loci: `61,894` (`6.19%`)
- Reads mapped to too many loci: `2,484` (discarded)

**Unmapped Reads**

- Unmapped Reads due to "too many mismatches": `0`
- Unmapped Reads due to "too short": `53,634` (`5.36%`)



### (2) Viewing BAM Files

#### ① Using `Samtools` on Linux

We can use `Samtools`, a widely used tool for *BAM* file operations, to quickly inspect alignment statistics, coverage depth, and specific region alignments.

**Installation:**

```bash
sudo apt update
sudo apt install samtools
```

`Samtools` can also be installed in a conda environment.

**Common Commands:**

1. **View BAM file contents:**

   ```bash
   samtools view 1M_SRR9336468_Aligned.sortedByCoord.out.bam | head
   ```

   ![Image](https://github.com/user-attachments/assets/5760a86f-0cf1-4de9-9a00-b9def2095dd2)

2. **Check alignment statistics:**

   ```bash
   samtools flagstat 1M_SRR9336468_Aligned.sortedByCoord.out.bam
   ```

   ![Image](https://github.com/user-attachments/assets/3498b1b7-0b95-4c90-bdf8-ec6339d6da11)

3. **Generate coverage statistics:**

   ```bash
   samtools depth 1M_SRR9336468_Aligned.sortedByCoord.out.bam > coverage.txt
   ```

4. **Generate BAM index file (`.bai`)**

   *BAM* index files speed up access to specific regions and improve processing efficiency. Tools like `IGV` require BAM index files for visualization.

   ```bash
   cd ~/STAR_results/
   samtools index 1M_SRR9336468_Aligned.sortedByCoord.out.bam
   ```

Additionally, we can visualize *BAM* files using the **IGV** software: [https://igv.org/doc/desktop/#DownloadPage/](https://igv.org/doc/desktop/#DownloadPage/)



#### ② Using Python Library `pysam`

**Installation**

Note: `pysam` does not run on Windows natively. We need to set up a *Python* environment on Linux using Conda. Note that the latest `pysam` version (0.22.1) does not yet support *Python* 3.13.

```bash
conda activate <env_name>
conda install pysam
```

![Image](https://github.com/user-attachments/assets/7ca3fa81-de9a-40dc-b36d-1fc22f55f7b5)

**Inspecting BAM Files with** `pysam`

**View File Contents:**

```python
import pysam

# Define BAM file path
bam_file_path = "../STAR_results/1M_SRR9336468_Aligned.sortedByCoord.out.bam"

# Open BAM file
with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
    print(f"Reading BAM file: {bam_file_path}")
    print(f"Reference genome: {bam_file.header['SQ']}")
    
    # Iterate through alignments
    for i, read in enumerate(bam_file.fetch()):
        print(f"Read {i + 1}:")
        print(f"  QNAME (Query Name): {read.query_name}")
        print(f"  FLAG: {read.flag}")
        print(f"  Reference Name: {read.reference_name}")
        print(f"  Start Position: {read.reference_start}")
        print(f"  Mapping Quality: {read.mapping_quality}")
        print(f"  CIGAR String: {read.cigarstring}")
        print(f"  Sequence: {read.query_sequence}")
        print(f"  Base Qualities: {read.query_qualities}")
        
        # Extract additional tags
        tags = dict(read.get_tags())
        nh = tags.get("NH")  # Number of hits
        hi = tags.get("HI")  # Hit index
        as_score = tags.get("AS")  # Alignment score
        nm = tags.get("NM")  # Number of mismatches

        print(f"  NH (Number of Hits): {nh}")
        print(f"  HI (Hit Index): {hi}")
        print(f"  AS (Alignment Score): {as_score}")
        print(f"  NM (Number of Mismatches): {nm}")
        print()
        
        # Limit output to first 5 reads
        if i >= 4:
            break
```

![Image](https://github.com/user-attachments/assets/d9516f9c-c9a8-44b4-b3fd-745472b3a74f)

**Check Alignment Statistics with** `pysam`

```python
# Initialize counters
total_reads = 0
primary_reads = 0
secondary_reads = 0
supplementary_reads = 0
duplicates = 0
mapped_reads = 0

# Open BAM file
with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
    for read in bam_file.fetch(until_eof=True):
        total_reads += 1
        if read.is_secondary:
            secondary_reads += 1
        elif read.is_supplementary:
            supplementary_reads += 1
        else:
            primary_reads += 1

        if read.is_duplicate:
            duplicates += 1

        if read.is_mapped:
            mapped_reads += 1

# Output results
print(f"{total_reads} total reads (QC-passed + QC-failed)")
print(f"{primary_reads} primary alignments")
print(f"{secondary_reads} secondary alignments")
print(f"{supplementary_reads} supplementary alignments")
print(f"{duplicates} duplicate reads")
print(f"{mapped_reads} mapped reads ({mapped_reads / total_reads * 100:.2f}%)")
```

![Image](https://github.com/user-attachments/assets/f8841d3c-65f8-4e46-8727-562a7f65acf3)
