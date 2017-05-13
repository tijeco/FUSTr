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
    identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",headerStr)
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
signature = "gene444::TR2|c0_g1_i135 len=243"
#signature = "rna0 gene=LOC102060483"
identifiers = re.search("c"+"(.*)"+"_g"+"(.*)"+"_i",signature)
if not Trinity_bool:
    print(signature[:identifiers.span()[1]].split("::")[1])
for i in range(len("ABC")):
    print(i)
#print(wordDict)

output = ["Families/family_26_dir/M2/statsfile.txt"]
identifiers = re.search("Families/"+"(.*)"+"_dir",output[0])
family=identifiers.groups()[0]
print(family)
# from pprint import pprint
# pprint(dir(identifiers))
# print(None)
#
# #print(identifiers.match())
# output = ["Families/family_2_dir/M2/statsfile.txt"]
# working_dir = output[0].split('/')[:-1][0] +'/'+output[0].split('/')[:-1][1]+'/'+output[0].split('/')[:-1][2]+'/'
# print("****************")
# print(working_dir)
# family = working_dir.strip("_dir/M2").strip("Families/")
#
# print(working_dir.strip("dir/M2"))
# print("apple".strip("ape"))
