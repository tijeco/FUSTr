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
SAMPLES2, = glob_wildcards("all.cds.combined_{sample}.fasta")


rule final:
    #input: "New.tmp"
    #input: expand("{sample}.pep.longestIsoform", sample=SAMPLES)
    input:expand("all.pep.combined_{sample2}.phy", sample2=SAMPLES2)
    #Aqinput:

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
        print("first ouput file",output.pep,"the following files")

        for i in input.pep_sequence:
            print(i)
        print("second ouput file",output.cds,"the following files")
        for i in input.cds_sequence:
            print(i)

        with open(output.pep, "w") as out:
            for i in input.pep_sequence:
                sample = i.split('.')[0]
                for line in open(i):
                    if ">" in line:
                        out.write(">"+sample+"_"+line.strip(">"))
                    else:
                        out.write(line)
        with open(output.cds, "w") as out:
            for i in input.cds_sequence:
                sample = i.split('.')[0]
                for line in open(i):
                    if ">" in line:
                        out.write(">"+sample+"_"+line.strip(">"))
                    else:
                        out.write(line)

                        #####Below shouldn't be necessary, but it might be ,
                        #### If it is then it will give one line sequences
"""
        with open(output.pep, "w") as out:
            for sample_file in input.pep_sequence:
                sample = sample_file.split('.')[0]
                sequence_iterator = fasta_iter(sample_file)
                for ff in sequence_iterator:
                    headerStr, seq = ff
                    out.write(">"+sample+"_"+headerStr)
                    out.write(seq)



"""



rule blastall:
    input:
        "all.pep.combined"
    output:
        "all.pep.combined.blastall.out"
    shell:
        " formatdb -i {input} -n {input}.seq.db;blastall -p blastp -d {input}.seq.db -i {input} -m 8 -o {output} -a 13"
rule mcl:
    input:
        "all.pep.combined.blastall.out"
    output:
        "all.pep.combined.mcl.dumpfile"
    shell:
        "mcxdeblast --m9 --line-mode=abc {input} -o {input}.abc;mcl {input}.abc --abc -I 2.0 -scheme 1 -o {output}"
rule mcl2tab:
    input:
        "all.pep.combined.mcl.dumpfile"
    output:
        "all.pep.combined_MCL.fnodes"
    run:
        number=1
        with open(input[0]) as f:
            with open(output[0], "w") as out:
                for line in f:
                    row = line.split()
                    for i in range(len(row)):

                        out.write(str(number)+"\t"+row[i]+"\n")
                    number+=1

rule sep_family_fasta:
    input:
        "all.pep.combined_MCL.fnodes"
    output:
        "TMP.file"
    shell:
        "cp all.cds.combined  all.cds.combined.fasta; silix-split -n 15 all.cds.combined.fasta {input} ; touch {output}"
#SAMPLES2, = glob_wildcards("MCL_CDS_FAM_15.members_dir/all.cds.combined_{sample}.fasta")

rule mafft_cds:
    input:
        "all.cds.combined_{sample2}.fasta"
    output:
        "all.cds.combined_{sample2}.aln"
    shell:
        "mafft --auto {input} > {output}"
rule mafft_pep:
    input :
        "all.pep.combined_{sample2}.fasta"
    output:
        "all.pep.combined_{sample2}.aln"
    shell:
        "mafft --auto {input} > {output}"

rule aln2phy:
    input:
        "all.pep.combined_{sample2}.aln"
    output:
        "all.pep.combined_{sample2}.phy"
    run:
        with open(output, "w") as out:

            sequence_iterator = fasta_iter(input)
            for ff in sequence_iterator:
                headerStr, seq = ff
                out.write(headerStr.strip('>').split()[0]+"\t")
                out.write(seq +"\n")
# rule mafft_tmpOneFile:
#     input:
#         expand("all.cds.combined_{sample}.aln", sample=SAMPLES)
#     output:
#         "New.tmp"
#     shell:
#         "touch {output}"
