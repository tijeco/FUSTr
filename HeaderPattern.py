import sys
import re

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-i" in sys.argv:
    fileName = getOptionValue("-i")
else:

    fileName = "test.fasta"
    #print("\nplease specify input file name using -i <file_name> \n")
    #sys.exit()
from itertools import groupby

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
Trinity_bool = False

#print(stringSplitter("d=';8\ouuiuuu.."))
sequence_iterator = fasta_iter(fileName)
fileLength = 0
wordDict = {}
subString = ""
pattern = ""
patternExists = True
for ff in sequence_iterator:
    fileLength+=1
    headerStr, seq = ff
    dentifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",headerStr)
    if identifiers !=None:
        Trinity_bool = True
    try:
        headerStr = headerStr.split()[0]+" " + headerStr.split()[1]
    except:
        headerStr = headerStr.split()[0]

    wordColumn = 1
    splitHeader = re.split(r'[`\ =~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', stringSplitter(headerStr))
    colNum = len(splitHeader)
    #print(splitHeader)
    try:
        usableColumns = min(colNum, usableColumns)
    except:
        usableColumns = colNum
    for i in range(usableColumns):
        #print(splitHeader[i])
        #print(i)
        try:
            wordDict[i][splitHeader[i]] = True
        except:
            wordDict[i] = {}
            wordDict[i][splitHeader[i]] = True
signature = ""

for i in range(usableColumns):
    #print(len(wordDict[i].keys()))
    if len(wordDict[i].keys()) == fileLength:
        #print("unique_id:", i)
        signature+="{unique_id}:"
    elif len(wordDict[i].keys()) == 1:
        for j in wordDict[i].keys():
            signature+=j + ":"
    else:
        signature+="{isoform_id}:"
    #print(i)
signature = signature[:-1]

print(signature)
signature = "TR2|c0_g1_i1 len=243"
signature = "rna0 gene=LOC102060483"
identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",signature)
print(identifiers==None)
if Trinity_bool:
    print("This is a trinity file!")
#print(wordDict)
