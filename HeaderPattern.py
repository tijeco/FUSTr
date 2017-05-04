import sys
import re

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-i" in sys.argv:
    fileName = getOptionValue("-i")
else:
    print("\nplease specify input file name using -i <file_name> \n")
    sys.exit()
from itertools import groupby

def fasta_iter(fasta_name):


    fh = open(fasta_name)


    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()[1:].strip()#Entire line, add .split[0] for just first column
        # print(header)


        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)
sequence_iterator = fasta_iter(fileName)
fileLength = 0
columnCountDict={}
wordDict = {}
newDict = {}
rowMembers = 1
subString = ""

RecentAlpha = False

for ff in sequence_iterator:
    fileLength+=1
    headerStr, seq = ff
    try:

        headerStr = headerStr.split()[0]+" " + headerStr.split()[1]
    except:
        headerStr = headerStr.split()[0]

    #print(headerStr.split()[0])
    wordColumn = 1
    print(re.split(r'[`\=~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', headerStr))
    for j in headerStr:

        try:
            if specialCharacterBool != (not j.isdigit() and not j.isalpha() and j!='-'):
                if wordColumn not in wordDict:
                    wordDict[wordColumn] = []
                    wordDict[wordColumn].append(subString)
                    newDict[wordColumn] = {}
                    newDict[wordColumn][subString] = True
                else:
                    newDict[wordColumn][subString] = True
                    if subString not in wordDict[wordColumn]:

                        wordDict[wordColumn].append(subString)


                #print wordColumn, subString
                wordColumn+=1
                #print subString
                subString = ""

            specialCharacterBool= (not j.isdigit() and not j.isalpha() and j!='-')
        except:
            specialCharacterBool= (not j.isdigit() and not j.isalpha() and j!='-')
        if specialCharacterBool:
            subString+=j
        else:
            subString += j
    print(wordColumn)
    if wordColumn not in wordDict:
        wordDict[wordColumn] = []
        wordDict[wordColumn].append(subString)
        newDict[wordColumn] = {}
        newDict[wordColumn][subString] = True
    else:
        newDict[wordColumn][subString] = True
        if subString not in wordDict[wordColumn]:
            wordDict[wordColumn].append(subString)


    subString= ""

print(wordDict)
print(newDict)

pattern= ""
numIsoformIDs = 0
for i in wordDict.keys():
    #print len(wordDict[i])
    if len(wordDict[i]) == 1:
        pattern+=wordDict[i][0]
    else:
        if len(wordDict[i]) == fileLength:

            pattern +="{unique_id}"
        else:
            pattern += "{isoform_id}"
            numIsoformIDs+=1
print(pattern)
patternDict = {}
for i in newDict.keys():
    print(len(newDict[i]))
    if len(newDict[i]) == 1:
        for j in newDict[i].keys():
            patternDict[i] = j
    elif len(newDict[i]) == fileLength:
        patternDict[i] = "{unique_id}"


print(patternDict)
#print(pattern)
#print(fileLength)
print(pattern)
