# FUSTr
Families Under Selection in Transcriptomes
# Introduction  
Fuster is a pipeline that clusters coding sequences from transcriptomes into protein families, and then analyzes those families for positive selection.

# Getting started

First you will need to install miniconda3, available  [here](https://conda.io/miniconda.html) so that you can set up the environment with all other dependencies. **Note:** make sure to choose the Python 3.6 version.

Alternatively, if you have pip3 installed, it is also possible to install miniconda using the following command.
```bash
pip3 install miniconda3
```
**Important:** The only dependency that is not avalaible in bioconda is SiLiX (v1.2.11), which needs to be installed on the system, and can be downloaded [here](http://lbbe.univ-lyon1.fr/Download,3009.html?lang=fr).

Once you have conda and SiLiX installed, clone into this repository.

```bash
git clone https://github.com/tijeco/FUSTr.git
```


Using conda, you can set up an environment that installs all dependencies in a local session

```bash
cd FUSTr
conda create --name FUSTr -c bioconda --file FUSTr.requirements.txt
source activate FUSTr
```




After that all you have to do is type the following command with the path to the directory containing all transcriptome assemblies (ending in .fasta) and the Snakefile will do the rest of the work!

```
snakemake -d <path_to_fastas>
```
Once that finishes you can end the conda session using

```bash
source deactivate  
```
