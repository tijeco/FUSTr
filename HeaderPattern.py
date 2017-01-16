
with open("sample.txt") as f:
    columnCountDict={}
    wordDict = {}
    rowMembers = 1
    for line in f:
        row = line.strip().split()
        columnCount = len(row)
        if len(row) not in columnCountDict:

            columnCountDict[len(row)] = 1
        else:
            columnCountDict[len(row)] += 1
        try:
            if len(columnCountDict)>rowMembers:
                print "columnCount has changed"
                rowMembers+=1
        except:
            None
        print columnCount
        #wordDict = {1:{1:"g",2:"e"},2:{1:"i",2:"d"}}
        numTypes = 0
        for i in range(columnCount):
            if i + 1 not in wordDict:
                wordDict[i+1] = {}
            for j in range(len(row[i])):
                # if j not in wordDict[i +1]:
                #     wordDict[i +1][j]=True
                if numTypes == 0:
                    numTypes = 1
                    # try:
                    #     int(row[i][j])
                    #     wordDict[numTypes] = {1:{1:{"g"}}}

                if wordDict[i+1]=={}:
                    wordDict[i+1][j]=row[i][j]
                try:
                    wordDict[i +1][j]
                except:
                    wordDict[i +1][j]=row[i][j]



print "this file has", rowMembers,"pattern/s per row"
for i in columnCountDict.keys():
    print columnCountDict[i],"row/s follow this pattern:",
    for j in range(i):
        print "{column"+str(j+1)+"}",
    print
print columnCountDict.keys()
print wordDict
print int("0")



from difflib import SequenceMatcher

string1 = "one_gene_one_"
string2 = "four_gene_four_"



match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
# if match.size == 0:
#     pattern = "{string}"
# else:
#     if match.a == 0 and match.b == 0:
#         if match.size == len(string1) and match.size == len(string2):
#             pattern = string1[match.a: match.a + match.size]
#         else:
#             pattern = string1[match.a: match.a + match.size]+" {string}"
#     else:
#         0


print(match)  # -> Match(a=0, b=15, size=9)
print(string1[match.a: match.a + match.size])  # -> apple pie
print(string2[match.b: match.b + match.size])  # -> apple pie

print string1[0:match.a]
print string1[match.a + match.size:len(string1)]
print string1.split(string1[match.a: match.a + match.size])
variableStrings1 = string1.split(string1[match.a: match.a + match.size])
variableStrings2 = string2.split(string2[match.b: match.b + match.size])
print variableStrings1
print variableStrings2
if len(variableStrings1) == 2:
    for i in range(len(variableStrings1)):
        match = SequenceMatcher(None, variableStrings1[i], variableStrings2[i]).find_longest_match(0, len(variableStrings1[i]), 0, len(variableStrings2[i]))
        print(string1[match.a: match.a + match.size])
