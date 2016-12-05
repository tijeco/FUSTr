# cds="ATGTTTGCTAAGATCGCTCTTCTGTGTGCTGCTATCGCAGTCGCTCAGTGCAACCCTGTCGTTCACGGACACAAAGCACTCGCTAGCACCGGAGTGAGTTCCAGGTCCCAGTCCCAAGATGGATACGGAAACTACGCTTTCGGTTATGACATCAAGGACGCTCTGGGTGCCACCAACTCC"
#
# original = "------------------------------------------------MFAKIALLCAAIAVAQCNPVVH--------------------GHKA----LASTGVSSRSQ-------SQD--GYGNYAFGYDIKDA---LGA-TNS------------------------------------------------------------------------------------------------------------"
# cut=[49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 138, 139, 142, 143, 144, 145, 146, 147, 148, 149]
#
# new = "FAKIALLCAAIAVAQCNLASTGVSSRSQYGNYAFGYDIKDALGTNS-----"
#
# CodonPos={}
# position=0
# codon=""
# number=1
# print len(cds)/3
# for i in cds:
#
#     codon +=i
#     #print i,position%3,codon
#     if position%3==2:
#         #print codon
#         #print codonTable[codon]
#         CodonPos[number]=codon
#         number+=1
#         #protein+=codonTable[codon]
#     position +=1
#
#     if position%3==0:
#         #print codon
#         codon=""
# print len(cut)
# print "****************"
# print CodonPos
# print "****************"
#
# aaPos = 1
# translated = ""
# for i in original:
#     if i !="-":
#         translated+=CodonPos[aaPos]
#     else:
#         translated+="-"
# print translated
# column=0
# trimmed=""
# aaPos=1
# prot=""
# for i in original:
#     if column  in cut:
#         if i =="-":
#             trimmed+="---"
#             prot+=i
#         else:
#             trimmed+=CodonPos[aaPos]
#             prot+=i
#             aaPos+=1
#     column+=1
# print column
# print trimmed
# print prot
# print len(prot)
# print len(trimmed)
# print len(trimmed)/3

original="-----------------GMRFALMEIKTCMAYVIAKFIIKKCPETRIPLEFKTSQGLLQPKSIILKMEIREDCPIRE--"

new="IAKFVIKRCPETKVPLEFNFGQGLLQPKEIILKVETREDSPIRE"

CodonPos= {1: 'GGA', 2: 'ATG', 3: 'AGG', 4: 'TTC', 5: 'GCC', 6: 'CTC', 7: 'ATG', 8: 'GAG', 9: 'ATC', 10: 'AAA', 11: 'ACG', 12: 'TGC', 13: 'ATG', 14: 'GCC', 15: 'TAC'
, 16: 'GTC', 17: 'ATT', 18: 'GCC', 19: 'AAG', 20: 'TTC', 21: 'ATC', 22: 'ATA', 23: 'AAG', 24: 'AAG', 25: 'TGT', 26: 'CCA', 27: 'GAA', 28: 'ACA', 29: 'CGG',
 30: 'ATT', 31: 'CCC', 32: 'TTG', 33: 'GAG', 34: 'TTT', 35: 'AAA', 36: 'ACT', 37: 'AGC', 38: 'CAA', 39: 'GGT', 40: 'TTG', 41: 'CTG', 42: 'CAA', 43: 'CCA',
44: 'AAA', 45: 'TCA', 46: 'ATA', 47: 'ATT', 48: 'CTA', 49: 'AAA', 50: 'ATG', 51: 'GAG', 52: 'ATC', 53: 'AGG', 54: 'GAA', 55: 'GAC', 56: 'TGT', 57: 'CCA',
58: 'ATT', 59: 'AGA', 60: 'GAA', 61: 'TGA'}

cut=[33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71
, 72, 73, 74, 75, 76]


aaPos=1
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

print prot
print trimmed

original="-----------------------------------------------------------------------------------MPE----------------QSNDYRVVVFGAGGVGKSSLVLRFVKGTFCEMYIPTVE-DTYRQVITCNRNI-YTLEITDTTGSHQFP-AMQRLNISKGHAFILVYAITSKQSLEELAPIFALIKEV--------------KGN---LEGIPMMLVGNKSDEEE------SREV-------------------------------------------------------------------------------------------------------"
original="-----------------------------------------------------------------------------------MPE----------------QSNDYRVVVFGAGGVGKSSLVLRFVKGTFCEMYIPTVE-DTYRQVITCNRNI-YTLEITDTTGSHQFP-AMQRLNISKGHAFILVYAITSKQSLEELAPIFALIKEV--------------KGN---LEGIPMMLVGNKSDEEE------SREV-------------------------------------------------------------------------------------------------------"


new = "FDLSKRSTFLSAMKWKKDVDACSSLPAILLGNKCDLQDRDV"
new = "YAITSKQSLEELAPIFALIKEGEGIPMMLVGNKSDEEEREV"
CodonPos={1: 'ATG', 2: 'CCC', 3: 'GAG', 4: 'CAA', 5: 'AGC', 6: 'AAT', 7: 'GAT', 8: 'TAC', 9: 'CGC', 10: 'GTG', 11: 'GTG', 12: 'GTG', 13: 'TTC', 14: 'GGA', 15: 'GCA'
, 16: 'GGG', 17: 'GGT', 18: 'GTG', 19: 'GGC', 20: 'AAG', 21: 'AGC', 22: 'TCC', 23: 'CTC', 24: 'GTG', 25: 'CTT', 26: 'CGC', 27: 'TTC', 28: 'GTC', 29: 'AAA',
 30: 'GGG', 31: 'ACC', 32: 'TTC', 33: 'TGC', 34: 'GAG', 35: 'ATG', 36: 'TAC', 37: 'ATT', 38: 'CCC', 39: 'ACG', 40: 'GTG', 41: 'GAG', 42: 'GAC', 43: 'ACC',
44: 'TAC', 45: 'AGG', 46: 'CAG', 47: 'GTG', 48: 'ATA', 49: 'ACG', 50: 'TGT', 51: 'AAC', 52: 'AGA', 53: 'AAT', 54: 'ATT', 55: 'TAC', 56: 'ACC', 57: 'CTC',
58: 'GAG', 59: 'ATC', 60: 'ACT', 61: 'GAT', 62: 'ACC', 63: 'ACA', 64: 'GGC', 65: 'AGC', 66: 'CAT', 67: 'CAG', 68: 'TTC', 69: 'CCA', 70: 'GCC', 71: 'ATG',
72: 'CAA', 73: 'AGG', 74: 'CTG', 75: 'AAC', 76: 'ATC', 77: 'AGT', 78: 'AAA', 79: 'GGA', 80: 'CAC', 81: 'GCC', 82: 'TTC', 83: 'ATA', 84: 'CTG', 85: 'GTG',
 86:'TAC', 87: 'GCC', 88: 'ATC', 89: 'ACC', 90: 'AGC', 91: 'AAA', 92: 'CAG', 93: 'AGC', 94: 'CTT', 95: 'GAA', 96: 'GAA', 97: 'CTG', 98: 'GCG', 99: 'CCA',
 100:'ATC', 101: 'TTC', 102: 'GCC', 103: 'CTC', 104: 'ATA', 105: 'AAG', 106: 'GAA', 107: 'GTG', 108: 'AAG', 109: 'GGC', 110: 'AAC', 111: 'CTG', 112: 'GAG',
 113: 'GGC', 114: 'ATC', 115: 'CCC', 116: 'ATG', 117: 'ATG', 118: 'CTG', 119: 'GTG', 120: 'GGC', 121: 'AAC', 122: 'AAG', 123: 'AGC', 124: 'GAC', 125: 'GAG',
126: 'GAG', 127: 'GAG', 128: 'AGC', 129: 'CGA', 130: 'GAG', 131: 'GTC'}

cut= [187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 224, 230, 231, 232, 233, 234, 235, 236, 237, 238,
 239, 240, 241, 242, 243, 244, 245, 253, 254, 255]


print original.count('-')
print len(original)-original.count('-')
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

print prot
print new
print trimmed
