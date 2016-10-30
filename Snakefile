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
SAMPLES, = glob_wildcards("{sample}.fasta")
rule final:
    input: expand("{sample}.headers.txt", sample=SAMPLES)

rule get_headers:
    input:
        "{sample}.fasta"
    output:
        "{sample}.headers.txt"
    shell:
        "cat {input} |grep ">" |sed -e 's/>//g' > {output}"
