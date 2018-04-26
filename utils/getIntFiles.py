import sys
import os.path
from itertools import groupby

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-d" in sys.argv:
    working_dir = getOptionValue("-d").strip("/")
else:

    print("\nplease specify input directory name using -d <directory_name> \n")
    sys.exit()


fusterID_file = working_dir+"/intermediate_files/fusterID.txt"
family_file = working_dir + "/final_results/famsUnderSelection.txt"
familyDir = working_dir + "/Families/"
selectDir = working_dir + "/final_results/famsUnderSelection_dir/"
os.makedirs( selectDir,  exist_ok=True)

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__()[1:].strip()#Entire line, add .split[0] for just first column

        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

id_dict = {}
# fams = {}
if os.path.exists(fusterID_file):
    if os.path.exists(family_file):
        #make id_dict
        with open(fusterID_file) as f:
            for line in f:
                row = line.strip().split()
                id_dict[row[0]] = row[1]
        with open(family_file) as f:
            line1 = True
            for line in f:
                if line1:
                    line1 = False
                else:
                    row = line.strip().split()
                    # fams[row[0]] = True
                    current_fam = row[0]

                    pep = familyDir + current_fam + ".fa"
                    cds = familyDir + current_fam + "_dir/" + current_fam + ".aln.codon"


                    new_pep = selectDir+ pep.split("/")[-1]
                    new_cds = selectDir+ cds.split("/")[-1]
                    print("opening:",cds)
                    print("writing to:",new_cds)

                    for current_file in [(new_pep,pep),(new_cds,cds)]:
                        with open(current_file[0],"w") as out:
                            sequence_iterator = fasta_iter(current_file[1])
                            for ff in sequence_iterator:
                                headerStr,seq = ff
                                new_header = id_dict[headerStr]
                                out.write(">"+new_header+"\n")
                                out.write(seq + "\n")







    else:
        print(family_file+" does not exist\nexiting")
        sys.exit()
else:
    print(fusterID_file +" does not exist\nexiting")
    sys.exit()
