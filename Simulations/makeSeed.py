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
    # length = random.randint(350, 500)
    length = random.randint(50, 1000)
    # branchlength = random.randint(1,15)
    branchlength = random.randint(1,50)
    # regime = random.choice(["pos","pur","con"])
    regime = "pos"
    filename = "seeds/seed"+regime+str(i)+"_len"+str(length)+"BL_"+str(branchlength)+".fa"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as out:
        out.write(makeSeed(length))
        # seed0_len100BL_10_True.tre
        # seed0_len100BL_10_True_alignment.FASTA
        # home/usr/hyphy/res/TemplateBatchFiles/FUBAR.bf
