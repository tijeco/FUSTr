# FUSTr
Families Under Selection in Transcriptomes
# Introduction  
Fuster is a pipeline that clusters coding sequences from transcriptomes into protein families, and then analyzes those families for positive selection.

# Getting started

The only software needed to run FUSTr is [Docker](https://www.docker.com/). FUSTr takes as input a directory containing transcriptome assemblies

**Important:** transcriptome assemblies need to be fasta files ending in .fasta, all in one directory.


Download FUSTr with the following command
```bash
git clone https://github.com/tijeco/FUSTr.git
```


With Docker installed correctly on your system issue the following command to initialize the Docker container

```
bash FUSTr/setup_docker.sh <directory_containing_fastas>
```

Once the docker container has been initialized you can enter it using the following command

```
docker run -it fustr /bin/bash
```

Now that you are in the docker container, your data is in /home/usr/data, to run FUSTr simply issue the following command
```
FUSTr -d ./data -t <number_of_threads>
```

The output will be in data/final_results/



You can use scp to transfer this to your local machine.

You do not have to make this Docker container over and over for new analysis with ```setup_docker.sh```. You can continue using the container for the analysis of additional datasets.


# Some notes about docker

In order to setup Docker on a new machine you will need root privileges for running commands or to create a group of users. This is not a problem if Docker is already properly installed on the system.

Also, the default container size for Docker is 10 GB, which was plenty to run the analysis for the manuscript (273,221 transcripts and 48,000 simulated transcripts). For larger datasets, this may not be enough space.

For the reasons above, in the event that users do not have root permissions to setup Docker on a new computer, or have a bewilderingly large dataset that would cause the 10GB Docker container to run out of space, below we have included instructions for installing FUSTr on the user's system.

# Installing FUSTr without Docker

These are the following dependencies of FUSTr that **must** be installed for FUSTr to properly function on a Linux 64 bit system.

1. [Miniconda3](https://conda.io/miniconda.html)
   * Be sure to choose **Python 3.6**
   * Will download ```Miniconda3-latest-Linux-x86_64.sh ``` for Linux 64 bit systems

   Install Miniconda3

   ```bash
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
   Choose to install to PATH
2.  [SiLiX](http://lbbe.univ-lyon1.fr/-SiLiX-?lang=en)

   * Download version **1.2.11**
   * Make sure Boost libraries are also installed (for Ubuntu issue the following commands, requires **root permissions**)
   ```bash
   sudo apt-get install libboost-dev
sudo apt-get install libboost-program-options-dev
   ```
   * Install SiLiX
```
tar zxvf silix-1.2.11.tar.gz
cd silix-1.2.11
./configure
make
make check
make install
```

3. [HYPHY](https://github.com/veg/hyphy)
   * Requires **root permissions**
   * First install ```cmake```
   ```bash
   apt-get install cmake
   ```
   * Install ```HYPHYMP```
   ```bash
   git clone https://github.com/veg/hyphy.git
   cd hyphy
   cmake .
   cmake .
   make HYPHYMP
   sudo make install
   ```
4. Add ```FUSTr``` to PATH
   * Add full path to FUSTr/bin to
    ```.bashrc``` file *i.e.*
    ```export PATH=/path/to/FUSTr/bin:$PATH```


   to run ```FUSTr``` simply issue the following command
   ```
   FUSTr -d directory_with_fastas -t <number_of_threads>
   ```



# What goes on under the hood

The file ```setup_docker.sh``` takes as input a directory that contains the transcriptome assemblies the user wishes to analyze. There may be any number of transcriptomes in this directory. The only reqirements are that they are
1. Uncompressed text files
2. Proper [FASTA format](https://en.wikipedia.org/wiki/FASTA_format)
3. End in .fasta
4. All contained in one directory

This directory is then added to a docker container (*located under /home/usr/data*) that installs all necessary third party dependencies, removing the need for users to install any of them on their actual system.

Once in the docker container (*automatically initiated at /home/usr*), ```FUSTr``` is installed to the system path.

Once ```FUSTr``` is executed using ```FUSTr -d ./data -t <number_of_threads>``` [Snakefile](http://snakemake.readthedocs.io/en/stable/) is executed to run 10 subroutines.

1. ```cleanFasta``` takes each  input fasta file in the ./data and cleans the text file for any spurious characters that commonly occur when transferring text files between different system architectures (such as ```^M```) that will break downstram analysis. The output for this file is found in **intermediate_files/{sample}.clean**

2. ```newHeaders``` takes as input the output from ```cleanFasta```. It further cleans the headers only keeping the first fields of text (some times headers have whole paragraphs of unnecessary information describing the sequence, so these are removed). Then these cleaned headers are analyzed to infer any patterns that may exist. The detected header patterns are placed in **headerPatterns.txt**, the output of the fasta files with new headers are placed in **intermediate_files/{sample}.new_headers**

3. ```orf``` takes as input the outputfrom ```newHeaders```. The program [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki) finds coding sequences from these transcripts, the output is placed in both **{sample}.new_headers.transdecoder.pep** and **{sample}.new_headers.transdecoder.cds**. Additional output from Transdecoder can be found in **{sample}.new_headers.transdecoder_dir**

4. ```longIsoform``` takes as input the files from ```orf``` and **headerPatterns.txt**
 as input to filter isoforms. It looks for genes that have multiple possible isoforms and only passes along the longest isoform for further analysis. The output can be found in **intermediate_files/{sample}.longestIsoform.pep** and **intermediate_files/{sample}.longestIsoform.cds**
5. ```blast``` takes the combined pep output from ```longIsoform``` with lighter unique identifiers as input using [DIAMOND](https://github.com/bbuchfink/diamond)
to run an all against all BLASTP search. The output can be found in **intermediate_files/all.pep.combined.blastall.out**

6. ```silix``` takes as input the output from ```blast``` and assigns proteins to putative gene families in the file **intermediate_files/all.pep.combined_r90_SLX.fnodes**

7. ```mafft``` generates multiple protein squence alignments for the families generated in ```silix``` into files **Families/family_{fam}.aln**
8. ```fasttree``` reconstructs phylogenies for the alignments generated by ```mafft``` into files **Families/family_{fam}.tree**

9. ```trimAln``` trims spurious columns from output of ```mafft``` in files **Families/family_{fam}.trimmed.aln**

10. ```hyphy``` takes trimmed alignments reverse translated to become codon alignments and classifies the selective regime of each site in each family. Output is placed in **Families/family_{fam}\_dir/family_{fam}.aln.codon.FUBAR.json**

# Investigating families of interest

Use the following command to parse only the json files listed in the **famsUnderSelection.txt** file and place them in a nice csv file per family. The columns will have the following information per codon position of the family alignment:

1. alpha Mean posterior synonymous substitution rate at a site
2. beta Mean posterior non-synonymous substitution rate at a site
3. beta-alpha Mean posterior beta-alpha
Prob[alpha>beta] Posterior probability of negative selection at a site
4. Prob[alpha<beta] Posterior probability of positive selection at a site
5. BayesFactor[alpha<beta] Empiricial Bayes Factor for positive selection at a site
6. PSRF Potential scale reduction factor - an MCMC mixing measure
7. Neff Estimated effective sample site for Prob [alpha<beta]

```bash
python FUSTr/utils/fubar_json.FUSTr.py -d <directory_with_fastas>
```
