from itertools import groupby
from itertools import (takewhile,repeat) #no longer needed 5.3.17
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML.chi2 import cdf_chi2
from Bio.Align.Applications import MafftCommandline #no longer needed 5.3.17
from io import StringIO #no longer needed 5.3.17
from Bio import AlignIO #no longer needed 5.3.17
from Bio import SeqIO #no longer needed 5.3.17
import sys #no longer needed 5.3.17
import re
from scipy import stats


def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__()[1:].strip()#Entire line, add .split[0] for just first column
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)
def stringSplitter(string):
    finalString = ""
    numSpecialChar = 0
    for i in string:
        specialCharacterBool= (not i.isdigit() and not i.isalpha() and i!='-')
        if specialCharacterBool:
            numSpecialChar+=1
        else:
            if numSpecialChar > 0:
                finalString+="_"
            finalString+=i
            numSpecialChar = 0
    return finalString


SAMPLES, = glob_wildcards("{sample}.fasta")
COUNTER = []


rule final:
    input:"finalStatsfile.txt"
    #input:dynamic("Families/family_{fam}_dir/family_{fam}.tree.fubar.csv")
    #input:dynamic("Families/family_{fam}_dir/M8a/tmp.txt")




rule cleanFasta:
    input:
        "{sample}.fasta"
    output:
        "{sample}.clean"
    run:
        Trinity_bool = False
        signature = ""
        sample = input[0].split('.')[0]
        sequence_iterator = fasta_iter(input[0])
        fileLength = 0
        columnCountDict={}
        wordDict = {}
        rowMembers = 1
        fileLength = 0
        newDict = {}
        subString = ""
        usableColumns = 0
        pattern = ""
        patternExists = True
        with open(output[0],"w") as out:
            for ff in sequence_iterator:

                headerStr, seq = ff
                trinity_identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",headerStr)
                if trinity_identifiers !=None:
                    Trinity_bool = True
                out.write(">"+headerStr+'\n')
                out.write(  seq +"\n")
                fileLength+=1
                splitHeader = re.split(r'[`\ =~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', headerStr)
                colNum = len(splitHeader)
                try:
                    usableColumns = min(colNum, usableColumns)
                except:
                    usableColumns = colNum
                for i in range(usableColumns):
                    try:
                        wordDict[i][splitHeader[i]] = True
                    except:
                        wordDict[i] = {}
                        wordDict[i][splitHeader[i]] = True
            for i in range(usableColumns):
                if len(wordDict[i].keys()) == fileLength:
                    signature+="{unique_id}:"
                elif len(wordDict[i].keys()) == 1:
                    for j in wordDict[i].keys():
                        signature+=j + ":"
                else:
                    signature+="{isoform_id}:"
            signature = signature[:-1]

            print(signature)


        with open("headerPatterns.txt","a") as out:
            #NOTE
                #1. if trinity, write file
                #2. if unique_id, try writing file
                    #*** a. if only one isoform_id, write file
                    #*** b. otherwise, just unique_id will be used in next rule
            if Trinity_bool:
                out.write(sample+"\t"+"TRINITY\n")
            elif "{unique_id}" in signature:
                out.write(sample+"\t"+signature+"\n")

rule newHeaders:
    input:
        "{sample}.clean"
    output:
        "{sample}.new_headers"
    run:
        try:
            patternDict = {}
            with open("headerPatterns.txt") as f:
                for line in f:
                    row = line.strip().split()
                    patternDict[row[0]] = row[1]
            with open(output[0],"w") as out:
                pattern = patternDict[input[0].split('.')[0]]

                if "{unique_id}" in pattern:
                    if pattern.count("{isoform_id}") == 1:
                        pattern_list = pattern.split(":")
                        for i in range(len(pattern_list)):
                            if pattern_list[i] == "{isoform_id}":
                                isoform_pos =  i
                                continue
                            elif pattern_list[i] == "{unique_id}":
                                unique_pos =  i
                                continue

                sequence_iterator = fasta_iter(input[0])
                for ff in sequence_iterator:

                    headerStr, seq = ff
                    if pattern == "TRINITY":
                        trinity_identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",headerStr)
                        new_header = stringSplitter(headerStr).split()[0]#[:trinity_identifiers.span()[1]+1]
                    if "{unique_id}" in pattern:
                        splitHeader = re.split(r'[`\ =~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', stringSplitter(headerStr))
                        if pattern.count("{isoform_id}") == 1:
                            new_header = splitHeader[unique_pos] + "___" + splitHeader[isoform_pos]
                        else:
                            new_header = splitHeader[unique_pos]
                    out.write( ">"+new_header+'\n')
                    out.write(seq+'\n')

        except:
            with open(output[0],"w") as out:

                sequence_iterator = fasta_iter(input[0])
                for ff in sequence_iterator:

                    headerStr, seq = ff
                    out.write( ">"+stringSplitter(headerStr).split()[0] +'\n')
                    out.write(seq+'\n')




rule transdecoderLongIsoforms:
    input:
        "{sample}.new_headers"
    output:
        "{sample}.new_headers.transdecoder_dir/longest_orfs.pep"
    shell:
            "TransDecoder.LongOrfs -t {input}  -m 30"

  
rule transdecoderPredict:
    input:
        fastaFile="{sample}.new_headers",LongOrfs="{sample}.new_headers.transdecoder_dir/longest_orfs.pep"
    output:
        "{sample}.new_headers.transdecoder.pep","{sample}.new_headers.transdecoder.cds"
    shell:
            "TransDecoder.Predict -t {input.fastaFile} --single_best_orf"









rule longestIsoformPep:
    input:
        "{sample}.new_headers.transdecoder.pep"
    output:
        "Temp/{sample}.longestIsoform.pep"
    run:

        with open(output[0], "w") as out:

            longIsoform={}

            sequence_iterator = fasta_iter(input[0])
            sample = input[0].split('.')[0]
            for ff in sequence_iterator:

                headerStr, seq = ff
                trinity_identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",headerStr)
                if trinity_identifiers != None:
                    # GeneID = headerStr[:trinity_identifiers.span()[1]].split("::")[1]
                    gene_header = headerStr.split("::")[1]
                    trinity_identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",gene_header)
                    GeneID = gene_header[:trinity_identifiers.span()[1]]
                else:
                    try:

                        GeneID=headerStr.split('___')[1].split('::')[0]
                    except:
                        reduced_header = stringSplitter(headerStr.split()[0].split("::")[0]+headerStr.split()[0].split("::")[1])
                        out.write('>'+sample+"_"+reduced_header+'\n')
                        out.write(seq + '\n')
                        continue
                if GeneID not in longIsoform:
                    longIsoform[GeneID] = [len(seq),headerStr,seq]
                else:
                    if longIsoform[GeneID][0] < len(seq):
                        longIsoform[GeneID] = [len(seq),headerStr,seq]
            for i in longIsoform.keys():
                out.write('>'+sample+'_'+i+'\n')
                # out.write('>'+sample+'_'+longIsoform[i][1].split("::")[0]+'\n')
                out.write(longIsoform[i][2]+'\n')
rule longestIsoformCDS:
    input:
        "{sample}.new_headers.transdecoder.cds"
    output:
        "Temp/{sample}.longestIsoform.cds"
    run:
        with open(output[0], "w") as out:

            longIsoform={}

            sequence_iterator = fasta_iter(input[0])
            sample = input[0].split('.')[0]
            for ff in sequence_iterator:

                headerStr, seq = ff
                trinity_identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",headerStr)
                if trinity_identifiers != None:
                    GeneID = headerStr[:trinity_identifiers.span()[1]].split("::")[1]

                try:

                    GeneID=headerStr.split('___')[1].split('::')[0]
                except:
                    reduced_header = stringSplitter(headerStr.split()[0].split("::")[0]+headerStr.split()[0].split("::")[1])
                    out.write('>'+sample+"_"+stringSplitter(headerStr.split()[0])+'\n')
                    out.write(seq + '\n')
                    continue
                if GeneID not in longIsoform:
                    longIsoform[GeneID] = [len(seq),headerStr,seq]
                else:
                    if longIsoform[GeneID][0] < len(seq):
                        longIsoform[GeneID] = [len(seq),headerStr,seq]
            for i in longIsoform.keys():

                out.write('>'+sample+'_'+longIsoform[i][1].split("::")[0]+'\n')
                out.write(longIsoform[i][2]+'\n')



rule combine_pep:
    input:
        expand("Temp/{sample}.longestIsoform.pep",sample=SAMPLES)
    output:
        "Temp/all.pep.combined","Temp/fusterID.txt"

    run:
        fusterID = 1
        with open(output[1],"w") as id_out:
            with open(output[0], "w") as out:
                for i in input:
                    sample = i.strip("Temp/").split('.')[0]
                    for line in open(i):
                        if ">" in line:
                            out.write(">fusterID_"+str(fusterID)+"\n")
                            id_out.write("fusterID_"+str(fusterID) + "\t"+line.strip(">"))
                            fusterID+=1
                        else:
                            out.write(line)
rule combine_cds:
    input:
        expand("Temp/{sample}.longestIsoform.cds",sample=SAMPLES)
    output:
        "Temp/all.cds.combined"

    run:
        fusterID = 1

        with open(output[0], "w") as out:
            for i in input:
                sample = i.strip("Temp/").split('.')[0]
                for line in open(i):
                    if ">" in line:
                        out.write(">fusterID_"+str(fusterID)+"\n")
                        fusterID+=1
                    else:
                        out.write(line)



rule blastall:
    input:
        "Temp/all.pep.combined"
    output:
        "Temp/all.pep.combined.blastall.out"
    shell:
        """
        makeblastdb -in {input} -out {input}.seq.db -dbtype prot
        blastp -db {input}.seq.db -query {input} -outfmt 6 -out {output} -num_threads 13 -evalue 1E-5
        """

rule silix:
    input:
        sequence_file="Temp/all.pep.combined",
        blast_file = "Temp/all.pep.combined.blastall.out"
    output:
        "Temp/all.pep.combined_r90_SLX.fnodes"
    shell:
        "silix -r 0.9 {input.sequence_file} {input.blast_file} > {output} || true"


rule node2families:
    input:
        node_file="Temp/all.pep.combined_r90_SLX.fnodes",
        sequence_file="Temp/all.pep.combined"
    output:
        dynamic("Families/family_{fam}.fa")
    run:
        famDict = {}
        seqDict={}
        print("opening",input.node_file)
        with open(input.node_file) as f:
            for line in f:
                row = line.split()
                if row[0] not in famDict:
                    famDict[row[0]]= [row[1]]
                else:
                    famDict[row[0]].append(row[1])
            #print(famDict)

        sequence_iterator = fasta_iter(input.sequence_file)
        for ff in sequence_iterator:
            headerStr, seq = ff
            seqDict[headerStr] = seq

        for i in famDict.keys():
            if len(famDict[i])>14:
                FileName = "Families/family_"+i+".fa"
                print(FileName)
                with open(FileName, "w") as out:
                    for j in famDict[i]:
                        out.write('>'+j+'\n')
                        out.write(seqDict[j]+'\n')
            else:
                FileName = "Families/family_"+i+".only_"+str(len(famDict[i]))+"_seqs.fas"
                with open(FileName, "w") as out:
                    for j in famDict[i]:
                        out.write('>'+j+'\n')
                        out.write(seqDict[j]+'\n')


rule mafft:
    input:
        "Families/family_{fam}.fa"
    output:
        "Families/family_{fam}.aln"
    shell:
        "mafft --auto --thread -1 {input} > {output}"

rule trimAln:
    input:
        "Families/family_{fam}.aln"
    output:
        trimmed_file="Families/family_{fam}.aln.trimmed",
        column_file="Families/family_{fam}.aln.trimmed.column_file"
    shell:
        "trimal -in {input} -out {output.trimmed_file} -gappyout -colnumbering > {output.column_file}"##

rule aln2phy:
    input:
        "Families/family_{fam}.aln",
        "Families/family_{fam}.aln.trimmed"
    output:
        "Families/family_{fam}.phy",
        "Families/family_{fam}.phy.trimmed"
    run:
        seq_length=0
        #print(output,"is output")
        #print(input,"is input")
        for currentFile in range(len(output)):
            print(output[currentFile],input[currentFile])

            with open(output[currentFile], "w") as out:


                sequence_iterator = fasta_iter(input[currentFile])
                first_line =True
                for ff in sequence_iterator:

                    headerStr, seq = ff
                    if first_line:
                        seq_length = len(seq)
                        num_lines = num_lines = sum(1 for line in open(input[currentFile]) if line[0]=='>')
                        out.write(str(num_lines)+" "+str(seq_length)+"\n")
                        first_line=False

                    seq_length = len(seq)
                    out.write(headerStr.strip('>')+"\t")
                    out.write(seq +"\n")

#print(longIsoform_CDS_combined)
rule phy2codon:
    input:
        untrimmed="Families/family_{fam}.phy",
        column_file="Families/family_{fam}.aln.trimmed.column_file",
        nucleotide="Temp/all.cds.combined"
    output:
        "Families/family_{fam}_dir/family_{fam}.codon.phylip",
        "Families/family_{fam}_dir/family_{fam}.aln.codon"

    run:
        cut = ""
        print(input.untrimmed)
        print(input.column_file)
        print(input.nucleotide)
        print(output)
        print("making cds dictionary")
        longIsoform_CDS_combined ={}
        sequence_iterator = fasta_iter(input.nucleotide)
        for ff in sequence_iterator:
            headerStr, seq = ff
            GeneID = headerStr
            if GeneID not in longIsoform_CDS_combined:
                    longIsoform_CDS_combined[GeneID] = seq
        #Open outout
        print(len(longIsoform_CDS_combined))
        with open(output[0], "w") as out:
            with open(output[1],"w") as out2:
                #Get  column cut file
                with open(input.column_file) as f:
                    for line in f:
                        cut  +=line.strip().split("#ColumnsMap")[1]
                    cut = cut.split(',')
                    cut = list(map(int, cut))
                #Get corresponding untrimmed Alignments, as original, line by line
                line1=True
                first_line=True
                with open(input.untrimmed) as f:
                    for line in f:
                        if line1:
                            line1=False
                            continue

                        row =line.strip().split()
                        original=row[1]#cds
                        header=row[0]

                        #NOTE, potetintal bug below, if exception then sequence isn't declared and it can't go forward, use continue probably
                        try:
                            sequence=longIsoform_CDS_combined[header]#original
                        except:
                            print(header,"not in dict")
                            continue
                        CodonPos={}
                        position=0
                        codon=""
                        number=1
                        for i in sequence:

                            codon +=i
                            if position%3==2:
                                CodonPos[number]=codon
                                number+=1
                            position +=1

                            if position%3==0:
                                codon=""
                        aaPos=0
                        firstAA=True
                        alnPos=0
                        prot=""
                        trimmed=""
                        for i in original:
                            if i!="-":
                                aaPos+=1

                            if alnPos in cut:
                                prot+=i
                                if i != "-":
                                    trimmed+=CodonPos[aaPos]
                                else:
                                    trimmed+="---"
                            alnPos+=1
                        num_lines = sum(1 for line in open(input.untrimmed) )
                        if first_line:
                            out.write(str(num_lines-1) + " " + str(len(trimmed)) + '\n')
                            first_line=False
                        out.write(header+'   '+trimmed+'\n')
                        out2.write(">"+header+"\n")
                        out2.write(trimmed+"\n")
rule FastTree:
    input:
        "Families/family_{fam}.aln"
    output:
        "Families/family_{fam}_dir/family_{fam}.tree"
    shell:
        "fasttree  -nosupport {input} > {output} || true"

rule copyTreeAln:
    input:
        tree_before="Families/family_{fam}_dir/family_{fam}.tree",
        aln_before="Families/family_{fam}_dir/family_{fam}.codon.phylip"
    output:
        treeM8a= "Families/family_{fam}_dir/M8a/family_{fam}.tree",
        treeM0="Families/family_{fam}_dir/M0/family_{fam}.tree",
        treeM1="Families/family_{fam}_dir/M1/family_{fam}.tree",
        treeM2="Families/family_{fam}_dir/M2/family_{fam}.tree",
        treeM3="Families/family_{fam}_dir/M3/family_{fam}.tree",
        treeM7="Families/family_{fam}_dir/M7/family_{fam}.tree",
        treeM8="Families/family_{fam}_dir/M8/family_{fam}.tree",
        alnM8a = "Families/family_{fam}_dir/M8a/family_{fam}.codon.phylip",
        alnM0="Families/family_{fam}_dir/M0/family_{fam}.codon.phylip",
        alnM1="Families/family_{fam}_dir/M1/family_{fam}.codon.phylip",
        alnM2="Families/family_{fam}_dir/M2/family_{fam}.codon.phylip",
        alnM3="Families/family_{fam}_dir/M3/family_{fam}.codon.phylip",
        alnM7="Families/family_{fam}_dir/M7/family_{fam}.codon.phylip",
        alnM8 = "Families/family_{fam}_dir/M8/family_{fam}.codon.phylip"
    shell:
        """cp {input.aln_before} {output.alnM8a}
        cp {input.aln_before} {output.alnM0}
        cp {input.aln_before} {output.alnM1}
        cp {input.aln_before} {output.alnM2}
        cp {input.aln_before} {output.alnM3}
        cp {input.aln_before} {output.alnM7}
        cp {input.aln_before} {output.alnM8}

        cp {input.tree_before} {output.treeM8a}
        cp {input.tree_before} {output.treeM0}
        cp {input.tree_before} {output.treeM1}
        cp {input.tree_before} {output.treeM2}
        cp {input.tree_before} {output.treeM3}
        cp {input.tree_before} {output.treeM7}
        cp {input.tree_before} {output.treeM8}
        """



###############################################################
"""

for PAML rule,

    have 2 outputs, M8 and M01237
    also an out put file with the following columns
        GeneFamily  ModelComparison ChiSq D.F p-value maybeFDR
    Also,, try to generate all BEB files that exist (a lot of times they don't)
        We'll use these in the final rule to plot familes with sites of strong selection


"""






#####################################################################3

rule M0:
    input:
        "Families/family_{fam}_dir/M0/family_{fam}.tree",
        "Families/family_{fam}_dir/M0/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M0/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M0_cml = codeml.Codeml()
        M0_cml.alignment = input[1]
        M0_cml.tree = input[0]
        M0_cml.out_file = output[0]
        M0_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M0_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M0_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M0_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M0_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M0_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M0_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M0_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M0_cml.set_options(model = 0)	         # models for codons:
        M0_cml.set_options(NSsites = [0])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M0_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M0_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M0_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M0_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M0_cml.set_options(fix_omega = 0)	     # 1: omega or omega_1 fixed, 0: estimate
        M0_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M0_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M0_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M0_cml.set_options(Small_Diff = .45e-6) # Default value.
        M0_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M0_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M0_results = M0_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M0_lnL = M0_results.get("NSsites").get(0).get("lnL")
            M0_np = len(M0_results.get("NSsites").get(0).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M0_lnL,M0_np,M0_cml.working_dir)
            with open(M0_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM0\t"+str(M0_np)+"\t"+str(M0_lnL)+"\n")
        except:
            with open(M0_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM0\tNA\tNA\n")



rule M1:
    input:
        "Families/family_{fam}_dir/M1/family_{fam}.tree",
        "Families/family_{fam}_dir/M1/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M1/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M1_cml = codeml.Codeml()
        M1_cml.alignment = input[1]
        M1_cml.tree = input[0]
        M1_cml.out_file = output[0]
        M1_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M1_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M1_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M1_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M1_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M1_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M1_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M1_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M1_cml.set_options(model = 0)	         # models for codons:
        M1_cml.set_options(NSsites = [1])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M1_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M1_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M1_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M1_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M1_cml.set_options(fix_omega = 0)	     # 1: omega or omega_1 fixed, 0: estimate
        M1_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M1_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M1_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M1_cml.set_options(Small_Diff = .45e-6) # Default value.
        M1_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M1_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M1_results = M1_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M1_lnL = M1_results.get("NSsites").get(1).get("lnL")
            M1_np = len(M1_results.get("NSsites").get(1).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M1_lnL,M1_np,M1_cml.working_dir)
            with open(M1_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM1\t"+str(M1_np)+"\t"+str(M1_lnL)+"\n")
        except:
            with open(M1_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM1\tNA\tNA\n")


rule M2:
    input:
        "Families/family_{fam}_dir/M2/family_{fam}.tree",
        "Families/family_{fam}_dir/M2/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M2/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M2_cml = codeml.Codeml()
        M2_cml.alignment = input[1]
        M2_cml.tree = input[0]
        M2_cml.out_file = output[0]
        M2_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M2_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M2_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M2_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M2_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M2_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M2_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M2_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M2_cml.set_options(model = 0)	         # models for codons:
        M2_cml.set_options(NSsites = [2])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M2_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M2_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M2_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M2_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M2_cml.set_options(fix_omega = 0)	     # 1: omega or omega_1 fixed, 0: estimate
        M2_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M2_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M2_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M2_cml.set_options(Small_Diff = .45e-6) # Default value.
        M2_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M2_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M2_results = M2_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M2_lnL = M2_results.get("NSsites").get(2).get("lnL")
            M2_np = len(M2_results.get("NSsites").get(2).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M2_lnL,M2_np,M2_cml.working_dir)
            with open(M2_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM2\t"+str(M2_np)+"\t"+str(M2_lnL)+"\n")
        except:
            with open(M2_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM2\tNA\tNA\n")

rule M3:
    input:
        "Families/family_{fam}_dir/M3/family_{fam}.tree",
        "Families/family_{fam}_dir/M3/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M3/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M3_cml = codeml.Codeml()
        M3_cml.alignment = input[1]
        M3_cml.tree = input[0]
        M3_cml.out_file = output[0]
        M3_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M3_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M3_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M3_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M3_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M3_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M3_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M3_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M3_cml.set_options(model = 0)	         # models for codons:
        M3_cml.set_options(NSsites = [3])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M3_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M3_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M3_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M3_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M3_cml.set_options(fix_omega = 0)	     # 1: omega or omega_1 fixed, 0: estimate
        M3_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M3_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M3_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M3_cml.set_options(Small_Diff = .45e-6) # Default value.
        M3_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M3_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M3_results = M3_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M3_lnL = M3_results.get("NSsites").get(3).get("lnL")
            M3_np = len(M3_results.get("NSsites").get(3).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M3_lnL,M3_np,M3_cml.working_dir)
            with open(M3_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM3\t"+str(M3_np)+"\t"+str(M3_lnL)+"\n")
        except:
            with open(M3_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM3\tNA\tNA\n")

rule M7:
    input:
        "Families/family_{fam}_dir/M7/family_{fam}.tree",
        "Families/family_{fam}_dir/M7/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M7/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M7_cml = codeml.Codeml()
        M7_cml.alignment = input[1]
        M7_cml.tree = input[0]
        M7_cml.out_file = output[0]
        M7_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M7_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M7_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M7_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M7_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M7_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M7_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M7_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M7_cml.set_options(model = 0)	         # models for codons:
        M7_cml.set_options(NSsites = [7])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M7_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M7_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M7_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M7_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M7_cml.set_options(fix_omega = 0)	     # 1: omega or omega_1 fixed, 0: estimate
        M7_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M7_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M7_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M7_cml.set_options(Small_Diff = .45e-6) # Default value.
        M7_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M7_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M7_results = M7_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M7_lnL = M7_results.get("NSsites").get(7).get("lnL")
            M7_np = len(M7_results.get("NSsites").get(7).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M7_lnL,M7_np,M7_cml.working_dir)
            with open(M7_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM7\t"+str(M7_np)+"\t"+str(M7_lnL)+"\n")
        except:
            with open(M7_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM7\tNA\tNA\n")

rule M8:
    input:
        "Families/family_{fam}_dir/M8/family_{fam}.tree",
        "Families/family_{fam}_dir/M8/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M8/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M8_cml = codeml.Codeml()
        M8_cml.alignment = input[1]
        M8_cml.tree = input[0]
        M8_cml.out_file = output[0]
        M8_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M8_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M8_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M8_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M8_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M8_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M8_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M8_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M8_cml.set_options(model = 0)	         # models for codons:
        M8_cml.set_options(NSsites = [8])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M8_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M8_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M8_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M8_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M8_cml.set_options(fix_omega = 0)	     # 1: omega or omega_1 fixed, 0: estimate
        M8_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M8_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M8_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M8_cml.set_options(Small_Diff = .45e-6) # Default value.
        M8_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M8_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M8_results = M8_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M8_lnL = M8_results.get("NSsites").get(8).get("lnL")
            M8_np = len(M8_results.get("NSsites").get(8).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M8_lnL,M8_np,M8_cml.working_dir)
            with open(M8_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM8\t"+str(M8_np)+"\t"+str(M8_lnL)+"\n")
        except:
            with open(M8_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM8\tNA\tNA\n")






rule M8a:
    input:
        "Families/family_{fam}_dir/M8a/family_{fam}.tree",
        "Families/family_{fam}_dir/M8a/family_{fam}.codon.phylip"
    output:
        "Families/family_{fam}_dir/M8a/family_{fam}.mcl"
    run:
        COUNTER.append(True)
        M8a_cml = codeml.Codeml()
        M8a_cml.alignment = input[1]
        M8a_cml.tree = input[0]
        M8a_cml.out_file = output[0]
        M8a_cml.working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'


        M8a_cml.set_options(noisy = 9)	         # 0,1,2,3,9: how much rubbish on the screen
        M8a_cml.set_options(verbose = 1)	     # 1: detailed output, 0: concise output
        M8a_cml.set_options(runmode = 0)	     # 0: user tree;  1: semi-automatic;  2: automatic
        M8a_cml.set_options(seqtype = 1)	     # 1:codons; 2:AAs; 3:codons-->AAs
        M8a_cml.set_options(CodonFreq = 2)	     # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        M8a_cml.set_options(clock = 0)	         # 0: no clock, unrooted tree, 1: clock, rooted tree
        M8a_cml.set_options(aaDist = 0)	         # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        M8a_cml.set_options(model = 0)	         # models for codons:
        M8a_cml.set_options(NSsites = [8])	     # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; Needs to be array
        M8a_cml.set_options(icode = 0)	         # 0:standard genetic code; 1:mammalian mt; 2-10:see below
        M8a_cml.set_options(Mgene = 0)	         # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
        M8a_cml.set_options(fix_kappa = 0)	     # 1: kappa fixed, 0: kappa to be estimated
        M8a_cml.set_options(kappa = 2)	         # initial or fixed kappa
        M8a_cml.set_options(fix_omega = 1)	     # 1: omega or omega_1 fixed, 0: estimate
        M8a_cml.set_options(omega = 1)	         # initial or fixed omega, for codons or codon-based AAs
        M8a_cml.set_options(getSE = 0)	         # 0: don't want them, 1: want S.E.s of estimates
        M8a_cml.set_options(RateAncestor = 0)	 # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
        M8a_cml.set_options(Small_Diff = .45e-6) # Default value.
        M8a_cml.set_options(cleandata = 0)	     # remove sites with ambiguity data (1:yes, 0:no)?
        M8a_cml.set_options(fix_blength = 0)	 # 0: ignore, -1: random, 1: initial, 2: fixed


        M8a_results=M8a_cml.run(verbose=True)
        identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
        family=identifiers.groups()[0]
        try:

            M8a_lnL = M8a_results.get("NSsites").get(8).get("lnL")
            M8a_np = len(M8a_results.get("NSsites").get(8).get("parameters").get("parameter list").split())
            print("@@@@@@@@@@@@@@@@@@@@")
            print(M8a_lnL,M8a_np,M8a_cml.working_dir+"statsfile.txt")
            with open(M8a_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM8a\t"+str(M8a_np)+"\t"+str(M8a_lnL)+"\n")
        except:
            with open(M8a_cml.working_dir+"statsfile.txt","w") as out:
                out.write(family+"\tM8a\tNA\tNA\n")




rule codemlModels:
    input:
        dynamic("Families/family_{fam}_dir/M0/family_{fam}.mcl"),
        dynamic("Families/family_{fam}_dir/M1/family_{fam}.mcl"),
        dynamic("Families/family_{fam}_dir/M2/family_{fam}.mcl"),
        dynamic("Families/family_{fam}_dir/M3/family_{fam}.mcl"),
        dynamic("Families/family_{fam}_dir/M7/family_{fam}.mcl"),
        dynamic("Families/family_{fam}_dir/M8/family_{fam}.mcl"),
        dynamic("Families/family_{fam}_dir/M8a/family_{fam}.mcl")
    output:
        "codmlLiklihood.txt"
    run:
        with open(output[0],"w") as out:
            for i in input:
                dir_1 = i.split('/')[0]
                dir_2 = i.split('/')[1]
                dir_3 = i.split('/')[2]

                for line in open(dir_1+"/"+dir_2+"/"+dir_3+"/"+"statsfile.txt"):
                    out.write(line)
rule revertFusterIDFams:
    input:
        fusterID = "Temp/fusterID.txt",
        family_file = "Temp/all.pep.combined_r90_SLX.fnodes"
    output:
        "FusterFamilies.txt"
    run:
        id_dict = {}
        with open(input.fusterID) as f:
            for line in f:
                row = line.strip().split()
                id_dict[row[0]] = row[1]
        with open(output[0],"w") as out:
            with open(input.family_file) as f:
                for line in f:
                    row = line.strip().split()
                    out.write(row[0]+"\t"+id_dict[row[1]]+"\n")

rule codmlStats:
    input:
        codemlFile="codmlLiklihood.txt",
        family_file="FusterFamilies.txt"
    output:
        "finalStatsfile.txt"
    run:
        ChiSq_dict = {}
        with open(output[0],"w") as out:
            with open(input.codemlFile) as f:
                for line in f:
                    row = line.strip().split()
                    #ChiSq_dict[row[0]] = True
                    try:
                        # print(row[0])
                        if row[0] not in ChiSq_dict:
                            ChiSq_dict[row[0]] = {}
                        try:
                            ChiSq_dict[row[0]][row[1]] = (float(row[2]),float(row[3]))

                        except:
                            ChiSq_dict[row[0]][row[1]] = (0,0)

                    except:
                        None

            # print(ChiSq_dict)

            for i in ChiSq_dict.keys():


                M3_M0_chiSq = 2*(ChiSq_dict[i]["M3"][1]-ChiSq_dict[i]["M0"][1])
                M3_M0_df = ChiSq_dict[i]["M3"][0]-ChiSq_dict[i]["M0"][0]
                M3_M0_pvalue = stats.chi2.sf(M3_M0_chiSq,M3_M0_df)
                # print(i,"M3_M0",M3_M0_pvalue)
                out.write(i+"\tM3_M0\t"+str(M3_M0_pvalue)+"\n")


                M2_M1_chiSq = 2*(ChiSq_dict[i]["M2"][1]-ChiSq_dict[i]["M1"][1])
                M2_M1_df = ChiSq_dict[i]["M2"][0]-ChiSq_dict[i]["M1"][0]
                M2_M1_pvalue = stats.chi2.sf(M2_M1_chiSq,M2_M1_df)
                # print(i,"M2_M1",M2_M1_pvalue)
                out.write(i+"\tM2_M1\t"+str(M2_M1_pvalue)+"\n")

                M8_M7_chiSq = 2*(ChiSq_dict[i]["M8"][1]-ChiSq_dict[i]["M7"][1])
                M8_M7_df = ChiSq_dict[i]["M8"][0]-ChiSq_dict[i]["M7"][0]
                M8_M7_pvalue = stats.chi2.sf(M8_M7_chiSq,M8_M7_df)
                # print(i,"M8_M7",M8_M7_pvalue)
                out.write(i+"\tM8_M7\t"+str(M8_M7_pvalue)+"\n")

                M8_M8a_chiSq = 2*(ChiSq_dict[i]["M8"][1]-ChiSq_dict[i]["M8a"][1])
                M8_M8a_df = ChiSq_dict[i]["M8"][0]-ChiSq_dict[i]["M8a"][0]
                M8_M8a_pvalue = stats.chi2.sf(M8_M8a_chiSq,M8_M8a_df)
                # print(i,"M8_M8a",M8_M8a_pvalue)
                out.write(i+"\tM8_M8a\t"+str(M8_M8a_pvalue)+"\n")
        #NOTE add here a print statement summarizing Families, number, avg size, singletons etc..





rule FUBAR:
    input:
        stat="finalStatsfile.txt",
        align="Families/family_{fam}_dir/family_{fam}.aln.codon",
        tree="Families/family_{fam}_dir/family_{fam}.tree"
    output:
        "Families/family_{fam}_dir/family_{fam}.tree.fubar.csv"
    shell:
        "(echo 1; echo 1;echo {input.align}; echo {input.tree}; echo 20;echo echo 5; echo 2000000; echo 1000000;echo 100;echo 0.5 )|HYPHYMP FUBAR.bf"
