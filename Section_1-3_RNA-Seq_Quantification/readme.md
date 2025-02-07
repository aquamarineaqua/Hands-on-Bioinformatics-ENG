**Chapter 1** RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification

**Section 1-3** RNA-Seq quantification

Adapted from *STAT 115 2021 Homework 1 Problem 4,5*

## Table of Contents
  - [0 Preparation: Creating a New Environment](#0-preparation-creating-a-new-environment)
  - [1 RNA-Seq Quantification](#1-rna-seq-quantification)
    - [(1) Installing RSEM](#1-installing-rsem)
    - [(2) STAR+RSEM](#2-starrsem)
      - [① STAR](#①-star)
      - [About the parameter: `--quantMode TranscriptomeSAM`](#about-the-parameter---quantmode-transcriptomesam)
      - [② Generating RSEM Index Directory](#②-generating-rsem-index-directory)
      - [③ RSEM Quantification](#③-rsem-quantification)
    - [(3) Salmon](#3-salmon)
      - [① Installing Salmon](#①-installing-salmon)
      - [② Generating Salmon Index Directory](#②-generating-salmon-index-directory)
      - [③ Running Salmon on Paired-End FASTQ Data](#③-running-salmon-on-paired-end-fastq-data)
    - [(4) Comparison of STAR+RSEM and Salmon Speed](#4-comparison-of-starrsem-and-salmon-speed)


## 0 Preparation: Creating a New Environment

In the previous session, we used `conda` to create a new environment named *star_env* and installed STAR, samtools, and other tools within it.

To facilitate the use of both **Linux-based tools (e.g., STAR, samtools)** and **programming in Python and R**, we will now create a new consolidated environment named `r_bio` (or whatever you like). This environment will include:

- **R and Python environments**, allowing coding in VS Code.
- **STAR, samtools, and other bioinformatics tools**, eliminating the need to switch between environments frequently.

The following Bash commands set up this environment with `r-base=4.4.1` and `python=3.10`.

> This step was mentioned in "Prerequisite 3." If you have already created such an environment, you can skip this section.

```
conda config --add channels defaults
conda config --add channels conda-forge  # Add conda-forge channel
conda config --add channels bioconda  # Add bioconda channel

conda create -n r_bio -c conda-forge r-base=4.4.1 python=3.10  # Create the r_bio environment with specified Python and R versions

conda activate r_bio  # Activate the environment

conda install -c conda-forge r-essentials  # Install essential R packages

conda install numpy pandas matplotlib openpyxl scipy sympy jupyter # Install commonly used Python libraries
conda install biopython
conda install pysam

conda install -c bioconda star  # Install STAR for sequence alignment
conda install -c bioconda samtools  # Install samtools for data processing
```



## 1 RNA-Seq Quantification

Transcript quantification plays a crucial role in RNA-seq data analysis. Numerous tools are available for quantifying expression at the transcript level. **`RSEM`** (*Bo Li et al., BMC Bioinformatics 2011*) is a software package that estimates the expression levels of *genes* and *isoforms* from single-end or paired-end RNA-Seq data. RSEM supports three different alignment tools: bowtie, bowtie2, or STAR.

**`Salmon`** (*Rob Patro et al., Nature Methods 2017*) is an ultra-fast, **alignment-free** method for expression quantification, which also corrects GC bias.

Here, we will use the previously selected high-quality alignment sample **1M_SRR9336471**, and perform quantification using both **STAR+RSEM** and **Salmon**.

### (1) Installing RSEM

```bash
conda activate r_bio
conda install -c bioconda rsem
```

### (2) STAR+RSEM

#### ① STAR

```bash
STAR --genomeDir ~/STAR_data/index \
     --readFilesIn ~/STAR_data/test_rnaseq/1M_SRR9336471_1.fastq ~/STAR_data/test_rnaseq/1M_SRR9336471_2.fastq \
     --runThreadN 4 \
     --quantMode TranscriptomeSAM \
     --outFileNamePrefix ~/STAR_results/1M_SRR9336471_ \
     --outSAMtype BAM SortedByCoordinate
```

#### About the parameter: `--quantMode TranscriptomeSAM`

By default, STAR generates a BAM file based on **genome**-aligned coordinates. Adding `--quantMode TranscriptomeSAM` instructs STAR to output an additional BAM file aligned to **transcriptome** coordinates (`Aligned.toTranscriptome.out.bam`). This file contains only **uniquely mapped reads** and maps reads to reference transcripts instead of genome coordinates. It serves as input for downstream quantification tools such as RSEM, Salmon, or Kallisto, enabling transcript-level quantification.

#### ② Generating RSEM Index Directory

Before running RSEM, we need to prepare the index using the `rsem-prepare-reference` command. This step processes the reference genome and gene annotation files (*FASTA* and *GTF*) into an RSEM-compatible format for subsequent expression quantification.

First, create a new folder to store the RSEM index:

```bash
mkdir -p ~/RSEM_data/index
```

Generate the RSEM index by specifying the GTF annotation file, genome FASTA file, and the output directory. The string `'rsem_index'` in the command represents the prefix for output files.

```bash
rsem-prepare-reference --gtf ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.gtf \
~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
~/RSEM_data/index/rsem_index
```

#### ③ RSEM Quantification

Create a directory to store the results:

```bash
mkdir -p ~/RSEM_results
```

Run the following command to perform quantification:

The `1M_SRR9336471_Aligned.toTranscriptome.out.bam` file from STAR is used as input for RSEM, while `~/RSEM_data/index/rsem_index` serves as the index. Since we are using paired-end RNA-Seq data, we must specify the `--paired-end` parameter. If omitted, RSEM defaults to single-end data processing. The results will be stored in `~/RSEM_results/RSEMOut`.

```bash
rsem-calculate-expression --no-bam-output \
--paired-end \
--time \
--bam \
-p 8 \
~/STAR_results/1M_SRR9336471_Aligned.toTranscriptome.out.bam \
~/RSEM_data/index/rsem_index \
~/RSEM_results/RSEMOut
```

The generated `RSEMOut.isoforms.results` file contains transcript-level information, including length, effective length, expected count, TPM, FPKM, and IsoPct.

![Image](https://github.com/user-attachments/assets/67be8167-039d-4dc2-af29-782f92421b0a)

### (3) Salmon

#### ① Installing Salmon

```bash
conda activate r_bio
conda install -c bioconda salmon
```

#### ② Generating Salmon Index Directory

**(1) Generating Transcriptome Sequence File**

To perform quantification using Salmon, we first need to generate the **transcriptome index** using the `salmon index` command. However, this requires a FASTA-format **transcriptome sequence file** first.

Note that the **genome sequence** and **annotation file** can be combined to generate the **transcriptome sequence**, as the annotation file defines the structure of each transcript in the genome (e.g., exon locations, order, and splicing patterns), while the genome sequence provides the actual DNA base sequences.

Therefore, we first extract the transcriptome sequence file from the genome sequence file (`Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa`) and the annotation file (`Saccharomyces_cerevisiae.R64-1-1.107.gtf`).

We use the **`gffread`** tool for this conversion. For more details, refer to the [official documentation](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread).

Install gffread:

```bash
conda install -c bioconda gffread
```

Generate the transcriptome sequence from the GTF file and genome sequence:

```bash
mkdir -p ~/Salmon_data/transcriptome_sequence
gffread ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.gtf \
        -g ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
        -w ~/Salmon_data/transcriptome_sequence/transcripts.fa
```

**(2) Generating the Transcriptome Index for Salmon**

Once the transcriptome sequence is obtained, we can create the transcriptome index.

Create a new directory to store the Salmon index:

```bash
mkdir -p ~/Salmon_data/index
```

Generate the Salmon index. The `-t` flag specifies the transcriptome sequence file, `-i` defines the output index directory, and `-k` sets the *k-mer* length (default is 31 and can be adjusted based on read length).

```bash
salmon index -t ~/Salmon_data/transcriptome_sequence/transcripts.fa \
-i ~/Salmon_data/index \
-k 31
```

> **Regarding Decoy Sequence Warning**
>
> During index generation, you might see the following warning:
>
> *The salmon index is being built without any decoy sequences. It is recommended that decoy sequences (either computed auxiliary decoy sequences or the genome of the organism) be provided during indexing.*
>
> This occurs because Salmon recommends including **decoy sequences**, typically the reference genome sequence or other non-transcribed regions. Decoy sequences help reduce mapping ambiguity and improve quantification accuracy, especially when handling non-specific reads. If decoy sequences are not included, transcript quantification might be affected by non-specific alignments, particularly in organisms with highly repetitive sequences or complex genomes.
>
> More details can be found in the [official documentation](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode).
>
> In this case, we proceed without using decoy sequences.

#### ③ Running Salmon on Paired-End FASTQ Data

```bash
salmon quant -i ~/Salmon_data/index \
             -l A \
             -1 ~/STAR_data/test_rnaseq/1M_SRR9336471_1.fastq \
             -2 ~/STAR_data/test_rnaseq/1M_SRR9336471_2.fastq \
             -p 4 \
             -o ~/Salmon_data/output
```

Parameter explanations:

- `-l A`: Auto-detect read orientation.
- `-1` and `-2`: Specify paired-end sequencing files R1 and R2 (for single-end sequencing, use `-r` instead).
- We continue using the **1M_SRR9336471** paired-end sequencing sample.
- `-p 4`: Number of threads used.
- `-o`: Output directory.

---

Once completed, the `quant.sf` file in the `output` directory contains transcript-level quantification results, including Length, EffectiveLength, TPM, and NumReads.

![Image](https://github.com/user-attachments/assets/7923cf07-bcc1-469e-81fa-69ab846831c5)

### (4) Comparison of STAR+RSEM and Salmon Speed

By checking the log files:

**Salmon** (`Salmon_data/output/logs/salmon_quant.log`):

```
start_time: [2025-01-28 15:24:59.905]
end_time: [2025-01-28 15:25:02.840]
```

**STAR** (`STAR_results/1M_SRR9336471_Log.final.out`):

```
                                 Started job on | Jan 27 14:39:06
                             Started mapping on | Jan 27 14:39:06
                                    Finished on | Jan 27 14:39:24
```

**RSEM** (`RSEM_results/RSEMOut.time`):

```
Aligning reads: 0 s.
Estimating expression levels: 14 s.
Calculating credibility intervals: 0 s.
```

- **Salmon**: Approximately **3 seconds**
- **STAR+RSEM**: **18 seconds** (STAR) + **14 seconds** (RSEM) = **32 seconds**
