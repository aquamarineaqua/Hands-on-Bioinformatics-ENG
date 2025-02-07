**Chapter 1** RNA-Seq Data Processing and Analysis: Alignment, Quality Control, and Quantification

**Section 1-2** RNA-Seq quality control

adapted from *STAT 115 2021 Homework 1 Problem 2,3*

## Table of Contents
  - [1 Performing QC with RSeQC](#1-performing-qc-with-rseqc)
    - [(1) Viewing Basic Statistics of a BAM File](#1-viewing-basic-statistics-of-a-bam-file)
    - [(2) RNA-seq Quality Control](#2-rna-seq-quality-control)
    - [(3) Gene Body Coverage Plot](#3-gene-body-coverage-plot)
  - [2. Plotting Gene Body Coverage Using Python or R](#2-plotting-gene-body-coverage-using-python-or-r)



After completing the **alignment**, the next step is to perform **quality control (QC)** on the alignment results.

In this section, we will use the `RSeQC` tool for QC, primarily evaluating **transcript integrity** and **genome coverage**.

## 1 Performing QC with RSeQC

The tool required is `RSeQC`. Official website: https://rseqc.sourceforge.net/.

> Note: RSeQC does not support Windows and must be installed on Linux.

Installation command:

```bash
pip install RSeQC
```

Or in a conda environment:

```bash
conda install RSeQC
```

![Image](https://github.com/user-attachments/assets/8b96bafd-fe3c-4c7b-b2d7-6686b95453d4)

Once installed, various command-line tools within `RSeQC` can be used. These tools are implemented as Python `.py` scripts.

### (1) Viewing Basic Statistics of a BAM File

In the previous section, we discussed using `samtools flagstat` and the `pysam` library to inspect *BAM* file alignment results. Now, we can also use `RSeQC` for *BAM* statistics:

```bash
bam_stat.py -i 1M_SRR9336468_Aligned.sortedByCoord.out.bam
```

![Image](https://github.com/user-attachments/assets/fcf554ec-e8d9-42e0-b4d4-607561ae78c6)

The `bam_stat.py` tool provides more granular statistical information, including:

- **mapq < mapq_cut (non-unique)**: Number of non-uniquely mapped reads with mapping quality (MAPQ) below the threshold.
- **mapq >= mapq_cut (unique)**: Number of uniquely mapped reads with MAPQ above the threshold.

### (2) RNA-seq Quality Control

**1. Obtain the Gene Model in `BED` Format**

To assess **transcript integrity**, we use `RSeQC` to calculate **TIN (Transcript Integrity Number)**.

This requires a gene model file in **`BED` format**. There are multiple ways to obtain a BED file. One method is using **Galaxy** ([Galaxy website](https://usegalaxy.org/)) to convert a genome **GTF file** into a **BED file**.

The figure below illustrates the conversion of *Saccharomyces_cerevisiae.R64-1-1.107.gtf* into a **BED12** format file.

![Image](https://github.com/user-attachments/assets/64494971-a95b-4b92-93ab-0cfec86a6856)

I have also uploaded this **BED file** (`Saccharomyces_cerevisiae.R64-1-1.107.bed`) to the repository.

Alternatively, you can convert the **GTF file** to a **BED file** manually. BED format is a simplified version of genome annotation, you can easily get resolution through Google.

**2. Generate the Corresponding BAM Index (`.bai`) Files**

Use `samtools` to generate `.bai` index files for BAM files:

```bash
cd ~/STAR_results/
samtools index 1M_SRR9336468_Aligned.sortedByCoord.out.bam
samtools index 1M_SRR9336471_Aligned.sortedByCoord.out.bam
samtools index 1M_SRR9336474_Aligned.sortedByCoord.out.bam
```

**3. Rename the Converted `.bed12` File to `.bed` and Run `tin.py`**

The command below runs `tin.py`. The `-i` argument can be a directory or a single BAM file.
Further information on the official documentation: [tin.py](https://rseqc.sourceforge.net/#tin-py)

```bash
tin.py -i ~/STAR_results -r ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.bed
```

> **TIN (Transcript Integrity Number)** and **MedTIN (Median Transcript Integrity Number)** are crucial metrics for assessing transcript integrity in RNA-seq data. These values help evaluate the quality and degradation level of RNA samples, particularly in studies analyzing RNA integrity.
>
> According to the [official RSeQC documentation](https://rseqc.sourceforge.net/#tin-py), **TIN** (Transcript Integrity Number) is analogous to the traditional **RIN** (RNA Integrity Number). While RIN is widely used at the **sample or transcriptome level** to assess RNA quality, it has limitations that TIN overcomes. TIN assigns a score to each transcript (range: **0 ≤ TIN ≤ 100**). **MedTIN** (the median TIN of all transcripts) serves as a metric for RNA integrity at the sample level.

**4. Output Files**

For each BAM file, two output files are generated:

- One ending in `summary.txt`
- One ending in `tin.xls`

![Image](https://github.com/user-attachments/assets/faff0914-cd43-4982-9a7f-49dbb4a43ecf)

---

Here are the **TIN scores** from the `summary.txt` files for three BAM alignment results:

| BAM File                                    | TIN (Mean) | TIN (Median) | TIN (Stdev) |
| ------------------------------------------- | ---------- | ------------ | ----------- |
| 1M_SRR9336468_Aligned.sortedByCoord.out.bam | 85.85      | 91.49        | 15.18       |
| 1M_SRR9336471_Aligned.sortedByCoord.out.bam | 88.09      | 92.94        | 13.58       |
| 1M_SRR9336474_Aligned.sortedByCoord.out.bam | 84.50      | 90.59        | 16.05       |



### (3) Gene Body Coverage Plot

The **Gene Body Coverage Plot** is a commonly used visualization tool in RNA-seq data analysis to evaluate sequencing uniformity. It illustrates the distribution of sequencing coverage across different positions of the **gene body**, helping assess data quality.

The `geneBody_coverage.py` tool in `RSeQC` can generate this plot. If at least three BAM files are provided as input, both a **line graph** and a **heatmap** will be generated.

For more details, refer to the official RSeQC documentation: [genebody-coverage.py](https://rseqc.sourceforge.net/#genebody-coverage-py)

```bash
cd ~/STAR_results/
geneBody_coverage.py -r ~/STAR_data/genome/Saccharomyces_cerevisiae.R64-1-1.107.bed -i 1M_SRR9336468_Aligned.sortedByCoord.out.bam,1M_SRR9336471_Aligned.sortedByCoord.out.bam,1M_SRR9336474_Aligned.sortedByCoord.out.bam -o output
```

![Image](https://github.com/user-attachments/assets/3ae850c4-254f-447c-9fad-a8d8a3ea8a06)

Upon completion, the following files will be generated:

![Image](https://github.com/user-attachments/assets/b2fbc6ea-d885-4a13-8c55-6c3ef78816b6)

Opening the *curves* and *heatmap* PDF files reveals the following:

![Image](https://github.com/user-attachments/assets/230798fb-226e-45a4-86c8-41f6c509bd26)

The alignment quality of the three samples is generally good, with *SRR9336474* and *SRR9336471* showing slightly better results than *SRR9336468*.

## 2. Plotting Gene Body Coverage Using Python or R

The `geneBody_coverage.py` script generates an output file `output.geneBodyCoverage.txt`, which contains the normalized read counts across gene body percentiles, as shown below:

![Image](https://github.com/user-attachments/assets/0153910e-0ed7-4cc0-947e-9958e41e4931)

We can use `Python` or `R` to visualize the **Gene Body Coverage** of each sample. The following scripts demonstrate this:

**Python**

```python
import pandas as pd
import matplotlib.pyplot as plt

with open('output.geneBodyCoverage.txt') as f:
    lines = f.readlines()

percentiles = {}
for i in range(1, len(lines)):
    # Skip i = 0, which is the header
    line_parts = lines[i].replace('\n', '').split('\t')
    file = line_parts[0]
    percs = [float(x) for x in line_parts[1:]]
    percentiles[file] = percs

df = pd.DataFrame(percentiles)
x = range(1, 101)

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
plt.suptitle('Gene Body Coverage', size=15, weight='bold')
for i in range(len(df.columns)):
    col = df.columns[i]
    print(col)
    axs[i].plot(x, df[col]/1000)
    axs[i].set_xlabel("Percentile of gene body (5' -> 3')")
    axs[i].set_ylabel("Read # (in thousands)")
    axs[i].set_title(
        'Library ' + col.split('_')[1]
    )

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
```

![Image](https://github.com/user-attachments/assets/af916265-3534-4811-bf9d-867d15f18534)

---

**R**

```R
# Set plot output size in Jupyter Notebook R environment
options(repr.plot.width = 15, repr.plot.height = 5)

# Read file
lines <- readLines("output.geneBodyCoverage.txt")

# Extract data
percentiles <- list()
for (i in 2:length(lines)) { # Start from line 2 (skip header)
  line_parts <- unlist(strsplit(lines[i], "\t"))
  file <- line_parts[1]
  percs <- as.numeric(line_parts[-1])
  percentiles[[file]] <- percs
}

# Convert to dataframe
df <- as.data.frame(percentiles)

# Define x-axis (1 to 100)
x <- 1:100

# Load ggplot2 for visualization
library(ggplot2)

# Set layout for multiple plots (1 row, 3 columns)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

# Generate plots
for (i in names(df)) {
  col <- i
  plot(
    x, df[[col]] / 1000,
    type = "l",
    col = "blue",
    lwd = 2,
    xlab = "Percentile of gene body (5' -> 3')",
    ylab = "Read # (in thousands)",
    main = paste("Library", strsplit(col, "_")[[1]][2])
  )
}
```

![Image](https://github.com/user-attachments/assets/8ef6ae0f-bc35-4d24-ac93-9a55bee96548)
