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



SAMPLES, = glob_wildcards("{sample}.pep")

rule final:
    input: "all.pep.combined"
    #input: expand("{sample}.pep.longestIsoform", sample=SAMPLES)


    #input: "all.combined.blastall.out"


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

rule subset_pep_and_cds:
    input:
        header="{sample}.longestIsoform.txt",
        sequence="{sample}.pep",
        cds_sequence = "{sample}.cds"
    output:
        pep="{sample}.pep.longestIsoform",
        cds="{sample}.cds.longestIsoform"
    shell:
        "cat {input.header} |awk '{{ print $1\"|\"$2 }}'|xargs faidx -f -d':' {input.sequence} >{output.pep}; cat {input.header} |awk '{{ print $1\"|\"$2 }}'|xargs faidx -f -d':' {input.cds_sequence} >{output.cds}"


rule combine_pep_and_cds:
    input:
        cds_sequence=expand("{sample}.cds.longestIsoform",sample=SAMPLES),
        pep_sequence=expand("{sample}.pep.longestIsoform",sample=SAMPLES)
    output:
        pep="all.pep.combined",
        cds="all.cds.combined"

    run:
        print("###################################")
        print("Output looks like this")
        for i in output:
            print(i)
        print(output)
        print("~~~~~~~~~~~~~~~~~~~~~~~~")
        print("INput looks like this")
        for i in input:
            print(i)
        print(input.pep_sequence)
        with open("all.pep.combined", "w") as out:
            out.write("Place holder text")
