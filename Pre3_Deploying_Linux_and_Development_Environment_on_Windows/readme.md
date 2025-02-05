Deploying Linux System and Development Environment on Windows

## 目录

  - [1. Installing Ubuntu via WSL](#1-installing-ubuntu-via-wsl)
  - [2 Installing Conda in WSL](#2-installing-conda-in-wsl)
    - [1. Download Miniconda for Linux](#1-download-miniconda-for-linux)
    - [2. Run the Installation Script](#2-run-the-installation-script)
    - [3. Follow the Installation Prompts](#3-follow-the-installation-prompts)
    - [4. Add `conda` to PATH](#4-add-conda-to-path)
  - [3 Managing Environments with Conda in Linux and Using Them in VS Code](#3-managing-environments-with-conda-in-linux-and-using-them-in-vs-code)
    - [Why Use Linux Conda Environments?](#why-use-linux-conda-environments)
    - [Configuring VS Code to Recognize WSL's Linux Conda Environments](#configuring-vs-code-to-recognize-wsls-linux-conda-environments)



When performing bioinformatics analysis, almost all commonly used tools run on Linux and operate via the command line for data processing and analysis. Additionally, certain Python libraries for bioinformatics are better supported on Linux.

If you are using a Windows system, you need a way to run **Linux** as well. Apart from traditional methods such as dual-boot installations and virtual machines, Windows now provides **WSL** (Windows Subsystem for Linux), allowing Windows users to simultaneously access both Windows and Linux.

This guide explains how to **deploy an Ubuntu system using WSL on Windows** (Ubuntu being the most commonly used Linux distribution) and set up the development environment.

## 1. Installing Ubuntu via WSL

Ensure your Windows system is up to date (latest versions of Windows 11 or Windows 10). Then, open **Command Prompt** or **PowerShell** and enter:

```powershell
wsl --install
```

This command will execute a series of operations, including downloading and installing the latest Ubuntu Linux distribution (as shown in the image below).

For more information, refer to the official tutorials:

- [https://learn.microsoft.com/en-us/windows/wsl/install](https://learn.microsoft.com/en-us/windows/wsl/install)
- [https://learn.microsoft.com/en-us/windows/wsl/setup/environment](https://learn.microsoft.com/en-us/windows/wsl/setup/environment)

![Image](https://github.com/user-attachments/assets/793d1101-62c5-4a7b-8579-7285ed009a6d)

Once Ubuntu is installed, restart Windows. Ubuntu will be added to the Start menu—click on the Ubuntu icon to launch it.

![Image](https://github.com/user-attachments/assets/c927211c-22e8-42b1-839b-bfaa219494ad)

After a short wait, set up your Linux username and password to successfully enter the Ubuntu system.

You will see that the latest version of Ubuntu (e.g., Ubuntu 24.04) has been installed. This interface is the Linux command line environment.

![Image](https://github.com/user-attachments/assets/b58665ee-d198-454b-848d-acb8497c43f7)

To update and upgrade software packages, run the following command upon first login:

```bash
sudo apt update && sudo apt upgrade
```

---

You can now access Ubuntu files from Windows. The default location of Ubuntu's file system in Windows is `\\wsl.localhost\Ubuntu`, which can be opened in File Explorer.

![Image](https://github.com/user-attachments/assets/68bf56db-d770-459d-b667-d36895a817d3)

Additionally, from Ubuntu's command-line interface, you can open the Windows File Explorer at the current Ubuntu directory using:

```
explorer.exe .
```

For more details on file system sharing, refer to the official documentation: [https://learn.microsoft.com/en-us/windows/wsl/filesystems#file-storage-and-performance-across-file-systems](https://learn.microsoft.com/en-us/windows/wsl/filesystems#file-storage-and-performance-across-file-systems)

---

We can also use **VS Code** for programming within WSL.

Navigate to the target directory in Ubuntu and enter:

```bash
code .
```

Allow access when prompted, and VS Code will open.

![Image](https://github.com/user-attachments/assets/9b6db745-ac1d-47e9-b593-51e60f353bec)

To create a new `.ipynb` notebook file, save it within the Linux directory. If the Linux directory does not appear in the save dialog, copy the path from the File Explorer and paste it into VS Code's file save dialog.

![Image](https://github.com/user-attachments/assets/55487394-3408-4954-a4cb-40526f803020)

![Image](https://github.com/user-attachments/assets/0be154c7-462f-48f1-b591-2b8ea4d1336f)

For more information on using VS Code with WSL, refer to: https://learn.microsoft.com/en-us/windows/wsl/tutorials/wsl-vscode

Now, you can code using *Python* Environment (conda) on Windows. The next section will explain how to manage environments in Linux using `conda` and configure *VS Code* to recognize and utilize conda environments within *WSL* for development.

## 2 Installing Conda in WSL

Since the Conda environment on *Windows* is separate from the one in *WSL*, and we need to use various bioinformatics tools on Linux, it is essential to install the Linux version of Miniconda or Anaconda in *WSL*. Here, we choose to install `Miniconda`.

With Conda, we can manage both *Python environments* and *APT packages* (such as bioinformatics tools like `STAR`, `samtools`, etc.) within Linux.

### 1. Download Miniconda for Linux

Run the following command in the Ubuntu terminal (the installation command can be found on Anaconda's official website):

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

### 2. Run the Installation Script

Once the installation file is downloaded to the current directory, execute the installation:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

### 3. Follow the Installation Prompts

- Read and agree to the *license agreement*.
- It is recommended to use the default installation path: `/home/<your-username>/miniconda3`.

### 4. Add `conda` to PATH

If the `conda` command is not recognized after installation, manually add *conda* to the environment variables.

Open the `.bashrc` file using a file manager.

![Image](https://github.com/user-attachments/assets/c0659f4b-7538-4743-b27d-40992a5db123)

Then append the following line:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
```

![Image](https://github.com/user-attachments/assets/0d4a9fe9-03f8-43ef-8750-0ff197127d71)

Save the file and then run the following command in the terminal:

```bash
source ~/.bashrc
```

![Image](https://github.com/user-attachments/assets/435fa7b3-c4d5-4c34-9ea6-01594b88b40c)

Afterward, check the installation by running:

```bash
conda --version
```

If the version number appears, the environment variable has been successfully added. Finally, restart the terminal.

## 3 Managing Environments with Conda in Linux and Using Them in VS Code

Now that we can manage environments using `conda` on Linux in *WSL*, we can also configure *VS Code* on *Windows* to recognize and use these Linux-based Conda environments for development.

### Why Use Linux Conda Environments?

As mentioned earlier, the Linux Conda environment allows management of both Linux applications (e.g., bioinformatics tools like `STAR`, `samtools`) and Python libraries. Some packages, such as `pysam` (used for parsing *BAM* files generated by *STAR*), are only supported on Linux.

For example, we can run the following commands in Linux to create a Conda environment that includes *Python*, *R*, and commonly used packages, along with necessary bioinformatics tools:

```bash
conda config --add channels defaults
conda config --add channels conda-forge  # Add channel
conda config --add channels bioconda  # Add channel

conda create -n r_bio -c conda-forge r-base=4.4.1 python=3.10  # Create new environment 'r_bio' with specific Python and R versions

conda activate r_bio  # Activate the environment

conda install -c conda-forge r-essentials  # Install common R packages

conda install numpy pandas matplotlib openpyxl scipy sympy jupyter  # Install common Python libraries
conda install biopython
conda install pysam

conda install -c bioconda star  # Install STAR for sequence alignment
conda install -c bioconda samtools  # Install samtools for data processing
```

Now, we can conduct bioinformatics analyses within the `r_bio` environment in Linux.

### Configuring VS Code to Recognize WSL's Linux Conda Environments

To enable *VS Code* to detect and utilize the Conda environments created in WSL:

1. Open VS Code on Windows.
2. Navigate to the Extensions Marketplace and search for `Remote Development`, then install it.

![Image](https://github.com/user-attachments/assets/afebf689-61d0-4287-9b44-16a89c3095cf)

Once installed, the remote resource manager in VS Code allows selecting *WSL* as the target environment. This enables using Linux's Conda environments directly, as VS Code will automatically detect them. See the following example:

![Image](https://github.com/user-attachments/assets/ef2f33c0-3d8b-4d13-89a9-1208a51790a7)

> **Note:** It is necessary to install additional extensions for the WSL environment, such as Jupyter, Python, and R. These can be installed directly from VS Code's extensions interface.

With these configurations, we can now develop efficiently in a Linux Conda environment within Windows using VS Code.
