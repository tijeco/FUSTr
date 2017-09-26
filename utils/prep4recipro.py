import sys
import os.path
from itertools import groupby

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__()[1:].strip()#Entire line, add .split[0] for just first column
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-d" in sys.argv:
    working_dir = getOptionValue("-d").strip("/")
else:

    print("\nplease specify input directory name using -d <directory_name> \n")
    sys.exit()
if "-i" in sys.argv:
    recipro_file = getOptionValue("-i").strip("/")
else:

    print("\nplease specify input directory name using -i <directory_name> \n")
    sys.exit()


pep_file = working_dir + "/intermediate_files/all.pep.combined"
fusterID_file  = working_dir + "/intermediate_files/fusterID.txt"

id_dict = {}
reciproDict = {}

pepDict = {}
if os.path.exists(fusterID_file):
    if os.path.exists(pep_file):
        #make id_dict
        with open(fusterID_file) as f:
            for line in f:
                row = line.strip().split()
                id_dict[row[1]] = row[0]
        with open(recipro_file) as f:
            for line in f:
                row = line.strip().split()
                num = row[0]
                header = row[1]
                if num not in reciproDict:
                    reciproDict[num] = {}
                    reciproDict[num][header] =True
                else:
                    reciproDict[num][header] = True
        with open(pep_file) as f:
            sequence_iterator = fasta_iter(pep_file)
            for ff in sequence_iterator:
                headerStr, seq = ff
                pepDict[headerStr] = seq
        for i in reciproDict.keys():
            with open("group."+i+".fa","w") as out:
                for j in reciproDict[i].keys():
                    # print(j,j in pepDict)
                    out.write(">"+j+"\n")
                    out.write(pepDict[id_dict[j]]+"\n")

                    #pepDict ---> {"fusterID":"seq"}
                    #reciproDict ---> {group:[header1,header2]}
                    #idDict  ---> {fusterID:header}
