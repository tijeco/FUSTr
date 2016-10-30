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
    input: expand("{sample}.longestIsoform.txt", sample=SAMPLES)
    #input: "all.combined.txt"

rule get_headers:
    input:
        "{sample}.fasta"
    output:
        "{sample}.headers.txt"
    shell:
        "cat {input} |grep '>' |sed -e 's/>//g' > {output}"

rule prep_headers:
    input:
        "{sample}.headers.txt"
    output:
        "{sample}.prepped_headers.txt"
    shell:
        "cat {input} |awk '{{print $1,$5}}'|cut -d':' -f1,3,8|awk -F':' '{{print $2,$3}}'|awk -F'|' '{{print $1,$2 }}' > {output}"




rule keep_longest_isoform:
    input:
        "{sample}.prepped_headers.txt"
    output:
        "{sample}.longestIsoform.txt"
    shell:
        "Rscript KeepLongestIsoform.R {sample}.prepped_headers.txt {sample}.longestIsoform.txt"

# rule combine_files:
#     input:
#         expand("{sample}.prepped_headers.txt", sample=SAMPLES)
#     output:
#         "all.combined.txt"
#     run:
#         with open(output[0], 'w') as out:
#             for i in input:
#                 sample = i.split('.')[0]
#                 for line in open(i):
#                     out.write(sample + ' ' + line)
