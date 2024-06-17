# Software to install on the AWS server <!-- omit from toc -->

>OBS! The software installation instructions provided here were only tested for a MacOS operating system with an Apple M2 chip.

## Table of contents <!-- omit from toc -->

- [BayPass](#baypass)
- [bcftools](#bcftools)
- [bedtools](#bedtools)
- [Delly2](#delly2)
- [Git](#git)
- [plink v1.9+](#plink-v19)
- [python3+ (library pandas)](#python3-library-pandas)
- [R](#r)
- [R packages](#r-packages)
- [snpEff](#snpeff)
- [Stacks](#stacks)
- [tabix](#tabix)
- [VCFtools](#vcftools)

## BayPass
>NOTE: BayPass is not working on the Mac machine after following the steps described below, not sure why...

Go to this [website](https://forgemia.inra.fr/mathieu.gautier/baypass_public/) and download the file. More details on installation can be found [here](https://forgemia.inra.fr/mathieu.gautier/baypass_public/-/blob/master/manual/BayPass_manual.pdf).

Install Xcode, https://stackoverflow.com/questions/71088796/how-to-install-gfortran-for-macos-monterey-version-12-2-1-m1 with `brew install gcc`.

Then this:
```bash
sudo xcodebuild -license accept # this did not work for me
```
And then checked in the terminal than gfortran is installed `gfortran --help`.

Note also that BayPass relies on OpenMP to implement multi-threading (usually installed on most Linux OS). For mac, use `openmp`:
```bash
brew install llvm
brew install libomp
```
For compilers to find libomp you may need to set:
```bash
# If you need to have llvm first in your PATH, run:
echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> /Users/apfuentesp/.bash_profile

# For compilers to find llvm and libomp you may need to set:
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
```
Install R dependencies (optional but recommended):
```R
install.packages(c("mvtnorm", "geigen", "poolfstat"))
```

To install the program:
```bash
# Uncompress the file
tar -xf baypass_public-master.tar.bz2

# clean up
rm baypass_public-master.tar.bz2

# go to the sources sub-directory
cd sources
make clean all FC=gfortran

# after compiling, one may run the command below to remove module-procedure and other output files that are not needed to run the executable
make clean

# give execution rights to the program if needed
#chmod +x g_baypass

# add the program to your path
echo 'export PATH="/Users/apfuentesp/Dropbox/Teaching/2024_physalia_adaptation_course/bioinfo-tools/baypass_public-master/sources/:$PATH"' >> ~/.bash_profile

source ~/.bash_profile

# test
g_baypass
```

## bcftools
Go to the [page](http://www.htslib.org/download/) and click on `bcftools-1.20` to download the file.

Using the Terminal, run these commands:
```bash
#  go to the directory where the file is located
cd /path/where/you/placed/the/file

# Uncompress the file
tar -xf bcftools-1.20.tar.bz2

# Clean up 
rm bcftools-1.20.tar.bz2

# Build the package from source
cd bcftools-1.20
./configure
make
make install  # if permission is denied, use `sudo make install` and type your mac password

# The executable programs will be installed to a bin subdirectory under your specified prefix, so you may wish to add this directory to your $PATH:
export PATH=/where/to/install/bin:$PATH    # for sh or bash users
# to know the path to your bin/ directory, in macos use `whereis bcftools`, and use that path
#export PATH=/usr/local/bin/bcftools:$PATH

# Check the package was correctly installed by calling the help page
bcftools -h
```

## bedtools
Follow the instructions given [here](https://github.com/arq5x/bedtools2/releases/tag/v2.31.1).


## Delly2
[Delly](https://github.com/dellytools/delly) is available as a statically linked binary, a singularity container (SIF file), a docker container or via Bioconda. You can also build Delly from source using a recursive clone and make.
```bash
git clone --recursive https://github.com/dellytools/delly.git

cd delly/

make all
```

## Git
Examine if Git is installed on your machine with `git --version`. If not, follow the instructions given [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) for the operating system available on your machine.


## plink v1.9+
Go to this [page](https://www.cog-genomics.org/plink/1.9/), and download the package for your operating system, preferably the stable build.

As this is a precompiled program, it is ready to go, no need for installation, just add it to the PATH.

```bash
# add the program to your path
echo 'export PATH="/Users/apfuentesp/Dropbox/Teaching/2024_physalia_adaptation_course/bioinfo-tools/plink_mac_20231211/:$PATH"' >> ~/.bash_profile

source ~/.bash_profile

# test if it works
plink
```

## python3+ (library pandas)
Examine if python3 is installed on your machine by typing in the Terminal `python3 --version`. If python3 is not installed, download it from [here](https://www.python.org/downloads/), and follow the instructions for the operating system of your machine.

To check if the library `pandas` is already installed, execute python by typing `python3`, then type `import pandas`. If not installed an error message will appear.

Check if pip3 is correctly installed with `pip3 --version`.

Upgrade your pip to avoid errors during installation. `pip3 install --upgrade pip`

Then add pip to the PATH by following these [instructions](https://graycode.ie/blog/how-to-add-python-pip-to-path/), in our case:
```bash
# Find out the path to python
python3 -m site --user-base
#/Users/apfuentesp/Library/Python/3.9

# for permanent solution, use that path here
echo 'export PATH="/Users/apfuentesp/Library/Python/3.9/bin:$PATH"' >> ~/.bash_profile
source ~/.bash_profile

# verify the modification
pip --version
```

To install the [pandas](https://pandas.pydata.org/docs/index.html#) library, use `pip install pandas`.

Check the library was installed correctly:
```bash
python3
import pandas  # if message appears, all good
```

## R
Follow the instructions given [here](https://cran.r-project.org/mirrors.html) for the operating system available on your machine.

## R packages
Open RStudio, and execute these commands:
```R
# for packages in CRAN
install.packages(c("ade4", "adegenet", "corrplot", "data.table", "dplyr", "geosphere", "ggVennDiagram", "gplots", "ggplot2", "magrittr", "reshape2", "RColorBrewer", "splitstackshape", "tibble", "tidyverse", "varhandle", "vcfR", "vegan"))

# for packages in Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!("devtools" %in% installed.packages())){install.packages(devtools)}

library(devtools)

BiocManager::install("edgeR")

BiocManager::install("LEA")

BiocManager::install("goseq")

devtools::install_github("petrelharp/local_pca/lostruct")

BiocManager::install("qvalue")

devtools::install_github("whitlock/OutFLANK")

BiocManager::install("StAMPP") # install via BioConductor

```
>Answer "yes" if asked to install additional libraries.

Check that the packages were installed correctly with:
```R
library(ade4)
library(adegenet)
library(corrplot)
library(data.table)
library(dplyr)
library(edgeR)
library(LEA)
library(lostruct)
library(geosphere)
library(ggplot2)
library(ggVennDiagram)
library(goseq)
library(gplots)
library(magrittr)
library(tibble)
library(tidyverse)
library(OutFLANK)
library(RColorBrewer)
library(reshape2)
library(splitstackshape)
library(StAMPP)
library(varhandle)
library(vcfR)
library(vegan)

```
If no error messages appear, then all good!

>If you get error messages, examine which additional packages might be missing/outdated, and need to be installed. For this, follow equivalent commands as shown before (directly from CRAN, Bioconductor, or a dedicated GitHub repository).

## snpEff

Follow the instructions given [here](https://pcingola.github.io/SnpEff/).

**Update java for at least version 11**
Instructions can be found [here](https://apple.stackexchange.com/questions/461260/how-to-fix-error-java-lang-unsupportedclassversionerror-on-macos). Then execute:
```bash
brew install openjdk@17
```
Now snpEff should run fine.

## Stacks
Download the package from this [website](https://catchenlab.life.illinois.edu/stacks/), and follow these instructions given [here](https://catchenlab.life.illinois.edu/stacks/manual/#install).

For Mac, it is necessary to make some adjustments before compilation, which are described [here](https://groups.google.com/g/stacks-users/c/QTPVpWeF5Vo):
```bash
# If one installed gcc using Homebrew with this command `brew install gcc`, the gcc files are stored in /opt/homebrew/bin/gcc-14 and /opt/homebrew/bin/g++-14, but stacks uses gcc and g++ (without the ..-14)
# Then, you need to ensure that 'gcc' and 'g++' refer to the versions installed by brew, by creating an alias for the commands without the ..-14 suffix:
cd /opt/homebrew/bin
ln -s gcc-14 gcc
ln -s g++-14 g++

# You can test if the right version of gcc is picked up by running
gcc --version
# it should display
#gcc (Homebrew GCC 14.1.0_1) 14.1.0

g++ --version
#g++ (Homebrew GCC 14.1.0_1) 14.1.0

# Once gcc and g++ were setup correctly, install stacks with
tar xfvz stacks-2.66.tar
# clean up
rm stacks-2.66.tar

cd stacks-2.66
./configure
make
sudo make install

# A default install will install files in the following way:
/usr/local/bin	# Stacks executables and Perl scripts.
```
The pipeline is now ready to run.

## tabix
>Perhaps no need to install it as it might come with bcftools

## VCFtools
Go to this [page](https://vcftools.github.io/downloads.html) and download the TAR version of the file. Then execute:
```bash
# uncompress file
tar -xvf vcftools-vcftools-v0.1.16-20-gd511f46.tar

# clean up
rm vcftools-vcftools-v0.1.16-20-gd511f46.tar

# for mac, install various missing packages using homebrew
#brew install autoconf
#brew install automake
#brew install pkg-config

# If you have trouble with zlib, for compilers to find zlib you may need to set:
#export LDFLAGS="-L/opt/homebrew/opt/zlib/lib"
#export CPPFLAGS="-I/opt/homebrew/opt/zlib/include"

# install
cd vcftools
./autogen.sh
./configure
make
make install # if permission is denied, use `sudo make install` and type your mac password

# test installation
man vcftools

```