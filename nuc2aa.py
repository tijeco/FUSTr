cds="ATGTTTGCTAAGATCGCTCTTCTGTGTGCTGCTATCGCAGTCGCTCAGTGCAACCCTGTCGTTCACGGACACAAAGCACTCGCTAGCACCGGAGTGAGTTCCAGGTCCCAGTCCCAAGATGGATACGGAAACTACGCTTTCGGTTATGACATCAAGGACGCTCTGGGTGCCACCAACTCC"

original = "------------------------------------------------MFAKIALLCAAIAVAQCNPVVH--------------------GHKA----LASTGVSSRSQ-------SQD--GYGNYAFGYDIKDA---LGA-TNS------------------------------------------------------------------------------------------------------------"
cut=[49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 138, 139, 142, 143, 144, 145, 146, 147, 148, 149]

new = "FAKIALLCAAIAVAQCNLASTGVSSRSQYGNYAFGYDIKDALGTNS-----"

CodonPos={}
position=0
codon=""
number=1
print len(cds)/3
for i in cds:

    codon +=i
    #print i,position%3,codon
    if position%3==2:
        #print codon
        #print codonTable[codon]
        CodonPos[number]=codon
        number+=1
        #protein+=codonTable[codon]
    position +=1

    if position%3==0:
        #print codon
        codon=""
print len(cut)
print "****************"
print CodonPos
print "****************"

aaPos = 1
translated = ""
for i in original:
    if i !="-":
        translated+=CodonPos[aaPos]
    else:
        translated+="-"
print translated
column=0
trimmed=""
aaPos=1
prot=""
for i in original:
    if column  in cut:
        if i =="-":
            trimmed+="---"
            prot+=i
        else:
            trimmed+=CodonPos[aaPos]
            prot+=i
            aaPos+=1
    column+=1
print column
print trimmed
print prot
print len(prot)
print len(trimmed)
print len(trimmed)/3
