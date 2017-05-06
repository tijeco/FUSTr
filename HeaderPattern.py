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
    try:

        headerStr = headerStr.split()[0]+" " + headerStr.split()[1]
    except:
        headerStr = headerStr.split()[0]

    wordColumn = 1
    splitHeader = re.split(r'[`\ =~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', stringSplitter(headerStr))
    #print(headerStr)
    #print(splitHeader)
    colNum = len(splitHeader)

    try:
        usableColumns = min(colNum, usableColumns)
    except:
        usableColumns = colNum
    # try:
    #     if colNum != usableColumns:
    #         #print("weirdness")
    #         usableColumns = min(colNum,usableColumns)
    # except:
    #     colNum = len(splitHeader)
    #
    #     usableColumns = colNum
    # colNum = len(splitHeader)

    #print(usableColumns,len(splitHeader))
    for i in range(usableColumns):
        #print(splitHeader[i])
        #print(i)
        try:
            wordDict[i][splitHeader[i]] = True
        except:
            wordDict[i] = {}
            wordDict[i][splitHeader[i]] = True


#print(wordDict)
signature = ""
for i in range(usableColumns):
    #print(len(wordDict[i].keys()))
    if len(wordDict[i].keys()) == fileLength:
        #print("unique_id:", i)
        signature+="{unique_id}:"
    if len(wordDict[i].keys()) == 1:
        for j in wordDict[i].keys():
            signature+=j + ":"
    else:
        signature+="{isoform_id}:"
signature = signature[:-1]

print(signature)

#     try:
#
#         if colNum != usableColumns:
#             print("Unable to detect isoforms")
#             patternExists = False
#             usableColumns = min(colNum,usableColumns)
#
#     except:
#         colNum = len(splitHeader)
#         usableColumns = colNum
#     colNum = len(splitHeader)
#     print(usableColumns,colNum)
#     if patternExists:
#         for j in headerStr:
#
#             try:
#                 if specialCharacterBool != (not j.isdigit() and not j.isalpha() and j!='-'):
#                     if wordColumn not in newDict:
#                         newDict[wordColumn] = {}
#                         newDict[wordColumn][subString] = True
#                     else:
#                         newDict[wordColumn][subString] = True
#                     wordColumn+=1
#                     subString = ""
#
#                 specialCharacterBool= (not j.isdigit() and not j.isalpha() and j!='-')
#             except:
#                 specialCharacterBool= (not j.isdigit() and not j.isalpha() and j!='-')
#             if specialCharacterBool:
#                 subString+=j
#             else:
#                 subString += j
#         if wordColumn not in newDict:
#             newDict[wordColumn] = {}
#             newDict[wordColumn][subString] = True
#         else:
#             newDict[wordColumn][subString] = True
#         subString= ""
#     else:
#         print("There is no pattern")
#
# pattern= ""
# numIsoformIDs = 0
# print(len(newDict.keys()))
# if patternExists:
#
#     for i in newDict.keys():
#         # print(len(newDict[i]))
#         if len(newDict[i]) == 1:
#             for j in newDict[i].keys():
#                 pattern+= j
#         elif len(newDict[i]) == fileLength:
#             pattern+= "{unique_id}"
#         else:
#             pattern+="{isoform_id}"
#     if pattern.count("{isoform_id}") >1:
#
#         if pattern == "{isoform_id}|{isoform_id}_{isoform_id}_{isoform_id} len={isoform_id}" or pattern == "TRINITY_{isoform_id}_{isoform_id}_{isoform_id}_{isoform_id} len={isoform_id}":
#             print("this is a trinity assembly")
#         else:
#             print("too many possible isoforms")
#     elif pattern.count("{isoform_id}") == 0:
#         print("No isoform indication detected")
#
# else:
#     print("this requires a pattern")
#
# print(pattern)
# print(usableColumns)
# print(newDict)
