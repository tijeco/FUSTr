from os.path import join
from itertools import groupby
def fasta_iter(fasta_name):


    fh = open(fasta_name)


    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()[1:].strip()
        # print(header)


        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)


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
    input: expand("{sample}.longestIsoform.fa", sample=SAMPLES)
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
        "Rscript KeepLongestIsoformID.R {input} {output}"


rule subset_fasta:
    input:
        header_subset="{sample}.longestIsoform.txt"
        #sequence_file = "{sample}.fasta"
    output:
        "{sample}.longestIsoform.fa"
    run:
        #fiter = fasta_iter({output})
        with open(output[0], 'w') as out:


            for i in input:
                sample = i.split('.')[0]
                fasta_file = sample+".fasta"
                fitter = fasta_iter(fasta_file)
                for ff in fiter:
                    headerStr, seq = ff
                for line in open(i):
                    ID = line.split()[0]+"|"+line.split()[1]
                    for ff in fitter:
                        headerStr,seq =ff
                        if ID in headerStr:
                            out.write(">"+headerStr)
                            out.write(seq)

        #for ff in fiter:

            #headerStr, seq = ff






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
