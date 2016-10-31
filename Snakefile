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

#SAMPLES, = glob_wildcards('{sample}.pep')
# print(SAMPLES)
# print(glob_wildcards('{sample}.pep'))
# print('{sample}.pep')
SAMPLES, = glob_wildcards("{sample}.pep")
rule final:
    #input: expand("{sample}.longestIsoform.cds", sample=SAMPLES)
    input: "all.combined.blastall.out"

rule get_headers:
    input:
        "{sample}.pep"
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


rule subset_pep:
    input:
        header="{sample}.longestIsoform.txt",
        sequence="{sample}.pep",
        cds_sequence = "{sample}.cds"
    output:
        pep="{sample}.longestIsoform.fa",
        cds="{sample}.longestIsoform.cds"
    shell:
        "cat {input.header} |awk '{{ print $1\"|\"$2 }}'|xargs faidx -f -d':' {input.sequence} >{output.pep}; cat {input.header} |awk '{{ print $1\"|\"$2 }}'|xargs faidx -f -d':' {input.sequence} >{output.cds}"
#
# rule subset_cds:
#     input:
#         header = "{sample}.longestIsoform.txt",
#         sequence = "{sample}.cds"
#     output:
#         "{sample}.longestIsoform.cds"



        # wanted = []
        # for i in input:
        #
        #     with open(i) as f:
        #         for line in f:
        #             line = line.strip()
        #             if line != "":
        #                 wanted.append(line)
        # with open(output[0], 'w') as out:
        #     for i in input:
        #         sample = i.split('.')[0]
        #         pep_file = sample+".pep"
        #         fitter = fasta_iter(pep_file)
        #         for ff in fitter:
        #             headerStr,seq =ff
        #             #with open(i) as f:
        #                 #for line in f:
        #             for line in wanted:
        #                 ID = line.split()[0]+"|"+line.split()[1]
        #                 if ID in headerStr:
        #                     out.write(">"+sample+"_"+headerStr+"\n")
        #                     out.write(seq+"\n")
# rule subset_cds:
#     input:
#         header_subset="{sample}.longestIsoform.txt"
#     output:
#         "{sample}.longestIsoform.cds"
#     run:
#         with open(output[0], 'w') as out:
#             for i in input:
#                 sample = i.split('.')[0]
#                 pep_file = sample+".cds"
#                 fitter = fasta_iter(pep_file)
#                 for ff in fitter:
#                     headerStr,seq =ff
#                     with open(i) as f:
#                         for line in f:
#                             ID = line.split()[0]+"|"+line.split()[1]
#                             if ID in headerStr:
#                                 out.write(">"+sample+"_"+headerStr+"\n")
#                                 out.write(seq+"\n")

# rule samtool_index:
#     input:
#         "{sample}.pep"
#     shell:
#         "samtools faidx {input}"
# rule rewrite_headers:
#     input:
#         "{sample}.pep"
#     output:
#         "{sample}.newHeader.pep"
#     shell:




rule combine_pep:
    input:
        expand("{sample}.longestIsoform.fa", sample=SAMPLES)
    output:
        "all.combined.pep"
    shell:
        "cat *longestIsoform.cds > all.combined.cds; cat *longestIsoform.fa > {output}"


rule split_pep:
    input:
        "all.combined.pep"
    output:
        "all.blastall.out"

    run:
        number = 1
        fitter = fasta_iter(input)
        for ff in variable:
            headerStr,seq =ff:


            fileName = number+"individual.fasta"

            with open(fileName. "wb") as out:
                out.write(headerStr+"\n")
                out.write(seq)
                number+=1


#
# rule blastall:
#     input:
#         "all.combined.pep"
#     output:
#         "all.combined.blastall.out"
#     shell:
#         "touch all.combined.blastall.out"
#         #"formatdb -i {input} -n all.combined.db; blastall -p blastp -d all.combined.db -i {input} -m 8 -o {output} "
#
#


#
# rule spectral_clust:
#     input:
#         "all.combined.blastall.out"
#     output:
#         "clusterx.out"
#     shell:
#         "../../../scps-0.9.8-Linux-amd64/bin/clusterx --param k_max=20 -t blast {input} -o {output}"
#
