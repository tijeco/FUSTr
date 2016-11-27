from os.path import join
from itertools import groupby
from itertools import (takewhile,repeat)
import sys
import datetime
today = datetime.date.today()
OrthoFinderDir = today.strftime('Results_%b%d')


# def getOptionValue(option):
#     optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
#     optionValue = sys.argv[optionPos + 1]
# if "--DIR" in sys.argv:
#     print(getOptionValue("--DIR"))
#     sys.exit()

def fasta_iter(fasta_name):


    fh = open(fasta_name)


    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()[1:].strip()
        # print(header)


        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)


def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen )


SAMPLES, = glob_wildcards("{sample}.pep")
SAMPLES2, = glob_wildcards("all.pep.combined_{sample}.fasta")
RESULTS, = glob_wildcards("sequenceDir/Resuls_{date}")
ORTHOGROUP, = glob_wildcards("Alignments/OG{orthogroup}.fa")

place4File = "sequenceDir/"+OrthoFinderDir+"/Alignments/OG{orthogroup}.out"
#print(expand("Alignments/OG{orthogroup}.phy",orthogroup=ORTHOGROUP))

rule final:
    #input:expand("Alignments/OG{orthogroup}.fa",orthogroup=ORTHOGROUP)

    input: "combined.txt"
    #input:expand("Alignments/OG{orthogroup}.fa",orthogroup=ORTHOGROUP)

    #input: expand("sequenceDir/"+OrthoFinderDir+"/Alignments/OG{orthogroup}.out", orthogroup=ORTHOGROUP)

    #input: expand("sequenceDir/{sample}.longestIsoform.pep.fasta", sample=SAMPLES)
    #input:expand("all.pep.combined_{sample2}.RAXML.out.tre", sample2=SAMPLES2)
    #Aqinput:

    #input: "all.pep.combined.blastall.out"


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
        #cds_sequence = "{sample}.cds"
    output:
        pep="{sample}.pep.longestIsoform",
        #cds="{sample}.cds.longestIsoform"
    shell:
        "cat {input.header} |awk '{{ print $1\"|\"$2 }}'|xargs faidx -f -d':' {input.sequence} >{output.pep}"#; cat {input.header} |awk '{{ print $1\"|\"$2 }}'|xargs faidx -f -d':' {input.cds_sequence} >{output.cds}"

rule longestIsoformDirectory:
    input:
        "{sample}.pep.longestIsoform"
    output:
        "sequenceDir/{sample}.longestIsoform.pep.fasta"
    shell:
        " cp {input} {output} "

rule listAlignments:
    input:
         "sequenceDir/"+OrthoFinderDir+"/Alignments/OG{orthogroup}.fa"
    output:
        "sequenceDir/"+OrthoFinderDir+"/Alignments/OG{orthogroup}.out"
    shell:
        "touch  {output}"

rule moveAlignments:
    input:
        expand("sequenceDir/{sample}.longestIsoform.pep.fasta",sample=SAMPLES)
    output:
        "Alignments/OG{orthogroup}.fa"
    shell:
        "cp sequenceDir/"+OrthoFinderDir+"/Alignments Alignments/"
        #"mkdir Alignments;cd sequenceDir/" +OrthoFinderDir+"/Alignments; for f in $(find . -maxdepth 1 -type f -exec sh -c 'test $( grep -c '>' {} | cut -f1 -d' ' ) -gt "+"14"+"' \; -print);do  cp  $f ../../../Alignments/$f;done"

rule aln2phy:
    input:
        "Alignments/OG{orthogroup}.fa"
    output:
        "Alignments/OG{orthogroup}.phy"
    run:
        seq_length=0
        print(output,"is output")
        print(input,"is input")
        with open(output[0], "w") as out:


            sequence_iterator = fasta_iter(input[0])
            first_line =True
            for ff in sequence_iterator:

                headerStr, seq = ff
                if first_line:
                    seq_length = len(seq)
                    num_lines = num_lines = sum(1 for line in open(input[0]) if line[0]=='>')
                    out.write(str(num_lines)+" "+str(seq_length)+"\n")
                    first_line=False

                seq_length = len(seq)
                out.write(headerStr.strip('>').split(':')[0]+"\t")
                out.write(seq +"\n")

rule mergePhys:
    input:
        "Alignments/OG{orthogroup}.phy"
    output:
        "combined.txt"
    run:
         with open(output[0], 'w') as out:
            for i in input:
                sample = i.split('.')[0]
                for line in open(i):
                    out.write(sample + ' ' + line)



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
        " makeblastdb -in {input} -out {input}.seq.db -dbtype prot ;blastp -db {input}.seq.db -query {input} -outfmt 6 -out {output} -num_threads 13 -evalue 1E-5"
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


rule raxml:
    input:
        "all.pep.combined_{sample2}.phy"
    output:
        "all.pep.combined_{sample2}.RAXML.out.tre"
    shell:
        "raxmlHPC-PTHREADS-AVX2 -p 18274 -m PROTGAMMAWAG -T 12 -# 1000 -s {input} -n {output}"

# rule mafft_tmpOneFile:
#     input:
#         expand("all.cds.combined_{sample}.aln", sample=SAMPLES)
#     output:
#         "New.tmp"
#     shell:
#         "touch {output}"
