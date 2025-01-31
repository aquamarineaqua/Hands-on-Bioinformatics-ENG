Managing R Environment with Anaconda and Writing R in Jupyter Notebook

## Table of Content
  - [Creating an R Environment in Anaconda](#creating-an-r-environment-in-anaconda)
    - [0 Issues with the Official Guide](#0-issues-with-the-official-guide)
    - [1 Create an R Environment](#1-create-an-r-environment)
    - [2 Install R Essentials bundle: r-essentials](#2-install-r-essentials-bundle-r-essentials)
    - [3 Write R in VS Code using Jupyter Plugin](#3-write-r-in-vs-code-using-jupyter-plugin)

For those who have been using Python for data analysis and machine learning, Anaconda is a familiar tool for managing different Python environments. Additionally, Jupyter Notebook (including the Jupyter plugin in VS Code) is a common choice for Python coding.

However, when switching to R, many users may find RStudio less intuitive for the following reasons:

1. The default R + RStudio setup installs and manages R packages in a global environment, lacking the environment management flexibility of Anaconda.
2. RStudio executes R scripts similarly to running a `.py` file, whereas Jupyter Notebook provides a more interactive, block-by-block execution and display interface.

This guide explains how to manage an R environment using Anaconda and write R scripts in Jupyter Notebook. You need to install Anaconda or Miniconda before following this guide.

This tutorial applies to both Windows and Linux.

## Creating an R Environment in Anaconda

### 0 Issues with the Official Guide

The official Anaconda guide (https://docs.anaconda.com/working-with-conda/packages/using-r-language/) suggests creating an R environment using:

```
conda create -n r_env r-essentials r-base
```

However, this approach has a limitation: it only installs R version 3.6, as the official Anaconda channel does not include newer R versions.

You can check available R versions in the terminal using:

```
conda search r-base
```

![Image](https://github.com/user-attachments/assets/c662bfc8-ec0f-4eb3-b155-405950f5781d)

From the output, the latest R versions are available in the `conda-forge` channel, whereas `pkg/r` does not provide the latest updates.

### 1 Create an R Environment

To install the latest R version, we must use the **conda-forge** channel.

Open Anaconda Prompt or a terminal and create an environment named `r_ds` (or choose whatever you like) with R version 4.4.1 from conda-forge:

```
conda create -n r_ds -c conda-forge r-base=4.4.1
```

**Note:** If you also need Python in this environment, specify the Python version explicitly. If not specified, Python will be installed automatically when `r-essentials` is installed, and the latest version will be used. For compatibility with bioinformatics tools on Linux, it is recommended to specify Python 3.10 explicitly:

```
conda create -n r_bio -c conda-forge r-base=4.4.1 python=3.10
```

Installation successful:

![Image](https://github.com/user-attachments/assets/c4eadaad-9862-4606-b315-b348ce757ac2)

### 2 Install R Essentials bundle: r-essentials

After creating the R environment, you can install R packages individually using:

```
conda install -c conda-forge r-package-name
```

Most R packages in Anaconda start with the prefix `r-`.

Alternatively, you can install the **R Essentials bundle** (`r-essentials`), which includes over 80 commonly used R packages such as IRKernel, dplyr, shiny, ggplot2, tidyr, caret, and nnet.

Reference for R packages in Anaconda: https://repo.anaconda.com/pkgs/r/. If a package exists in this repository, it can be installed using: `conda install -c conda-forge r-package-name` .

To install `r-essentials`, first activate the environment:

```
conda activate r_ds
conda install -c conda-forge r-essentials
```

Once `r-essentials` is installed, Jupyter Notebook will recognize the **R kernel** (Jupyter installed as well), allowing R script execution within Jupyter.

![Image](https://github.com/user-attachments/assets/04336901-1e36-4a5b-9f18-8089c1d54be0)

With this setup, now Jupyter Notebook can run R scripts. Open Jupyter Notebook by executing:

```
jupyter notebook
```

![Image](https://github.com/user-attachments/assets/afe63b7a-99a9-417e-a8f9-dd177f646985)

![Image](https://github.com/user-attachments/assets/8c0b4512-e93f-4dba-959f-72e184e4bce0)

### 3 Write R in VS Code using Jupyter Plugin

If you love the Jupyter plugin in VS Code for Python, you may also want to write R in the same way.

Since we have already installed the `r_ds` environment in Anaconda, follow these steps in VS Code:

1. Open VS Code and install the R extension (ensure the Jupyter extension is installed first).

   ![Image](https://github.com/user-attachments/assets/46fdceed-ee95-4a4f-b06d-37450a1acc62)

2. Create a `.ipynb` file and click the environment selector in the top-right corner.

   ![Image](https://github.com/user-attachments/assets/2b3d4da8-d61f-4831-81a4-9e5347f42880)

3. Choose "Jupyter Kernel" and select the `r_ds` environment.

   ![Image](https://github.com/user-attachments/assets/31b8ec68-501e-46f8-8b6e-b9b896a5b309)

Now, you can enjoy R in VS Code on Jupyter Notebook!

---

Here’s my example of writing R in VS Code's Jupyter Notebook. Cool!

![Image](https://github.com/user-attachments/assets/692a9b24-da79-43b5-b146-613b3111dc7a)
