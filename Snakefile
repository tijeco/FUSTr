from os.path import join

# Globals ---------------------------------------------------------------------

# Full path to a FASTA file.
#GENOME = 'genome.fa'

# Full path to a folder that holds all of your FASTQ files.
#FASTQ_DIR = 'fastq/'

#SAMPLES, = glob_wildcards('{sample}.fasta')
# print(SAMPLES)
# print(glob_wildcards('{sample}.fasta'))
# print('{sample}.fasta')
rule all:
    input: "test.txt"

rule downsampling:
    input:
        "{bname}.fasta"
    output:
        "{bname}.txt"
    shell:
        "cat {input}  > {output}"
