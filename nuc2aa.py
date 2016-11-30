nuc ="ATGAACTTCTCTACAATCATCATACTACTGATTATTGGCATTTCTTACGCTACTGGACTGTCAACTAATCAAGGATCAACCAGAAGTAAAAAAGGAGCTCTTATGGAAATTGAAAGAAGGTGTATTGGAAAAGAGCAGGAATGCACTGATAACAAAGGAGGGTGCTGTGGTGACATGAAATGCTTATGTTACAAGAATGTCATGGATGAAAACGAGAAAGGATGCTGGTGTGGATCTGGAACCTACTACTTCATTGAAAAACCTTAA"
print len(nuc)/3
position = 0
codon=""
codonTable = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"C", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
protein =""
for i in nuc:
    codon +=i
    print i,position%3,codon
    if position%3==2:
        print codon
        print codonTable[codon] 
        protein+=codonTable[codon]
    position +=1

    if position%3==0:
        #print codon
        codon=""
print protein
print "MNFSTIIILLIIGISYATGLSTNQGSTRSKKGALMEIERRCIGKEQECTDNKGGCCGDMKCLCYKNVMDENEKGCWCGSGTYYFIEKP*"

original = "------------------MNFSTIIILLIIGISYATGLSTNQGSTRSKKGALMEI-ERRCIGKEQECTDNKGGCCG-----DMKCLCYKN-----VMDENEKGCWCGSGTYYFIE-KP"
print len(original),original
trimmed= "MNFSTIIILLIIGISYATGLSTNQGSTRSKKMEI-ERRCIGKEQECTDNKGGCCGDMKCLCYKNDENEKGCWCGSGTYY"
