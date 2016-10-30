# FUSTr
Families Under Selection in Transcriptomes
#Introduction
Fuster is a pipeline that clusters coding sequences from transcriptomes into protein families, and then analyzes those families for positive selection.

#Getting started

First you will need [pip3](http://stackoverflow.com/questions/6587507/how-to-install-pip-with-python-3/6587528#6587528) installed.
Using pip3 you can install miniconda3 in order to get conda so that you can set up the environment with all other dependencies.
"""bash
pip3 install miniconda3
"""
This program works on by analyzing transcriptome files that have been processed with transdecoder.
It can work on one or as many transcriptomes as you desire.


Once you have conda installed,change in to the directory that has your transcriptomes and clone into this repository.

"""bash
git clone https://github.com/tijeco/Fuster.gitThe best practice for running this program would be to have all of your transcriptome fasta files in one directory.
"""

The best practice for running this program would be to move all of your transcriptome fasta files into the Fuster directory.

Using conda, you can set up an environment that installs all dependencies in a local session

"""bash
cd Fuster
conda create --name Fuster -c bioconda --file Fuster.requirements.txt
source activate Fuster
"""

Now all you have to do is type 'snakemake'
