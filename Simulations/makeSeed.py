import sys
import random
import os
def makeCodon():
    codon = ""
    AA = "ATCG"
    # print(len(AA),print(AA[0],AA[19]))
    # print( random.randint(0,len(AA)))
    codon += AA[random.randint(0,len(AA)-1)]
    codon += AA[random.randint(0,len(AA)-1)]
    codon += AA[random.randint(0,len(AA)-1)]
    return codon
# print(makeCodon())
def makeSeed(numAA):
    seed ="ATG"
    for i in range(numAA):
        codon = makeCodon()
        if codon not in ["TAA","TAG","TGA"]:
            seed+=codon
    seed+="TAA"
    return seed

# seedNum = int(sys.argv[1])

for i in range(int(sys.argv[1])):
    length = random.randint(50, 150)
    branchlength = random.randint(1,15)
    regime = random.choice(["pos","pur","con"])
    filename = "seeds/seed"+str(i)+"_len"+str(length)+"BL_"+str(branchlength)+"."+regime+".fa"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as out:
        out.write(makeSeed(length))
        # seed0_len100BL_10_True.tre
        # seed0_len100BL_10_True_alignment.FASTA
        # home/usr/hyphy/res/TemplateBatchFiles/FUBAR.bf
