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




family_file = working_dir + "/final_results/FUSTr.fams"

species_dict = {}

num_fams = 0
with open(family_file) as f:
        for line in f:
                row = line.strip().split()
                num = int(row[0])
                species = row[1].split("_")[0]
                if species not in species_dict:
                        species_dict[species] = {}
                species_dict[species][num] = True
                num_fams+=1
line2print = ""
print(str(len(species_dict))+"       "+str(num_fams))
for i in species_dict.keys():
        line2print = ""
        for j in range(1,int(num_fams)+1):
                if j in species_dict[i]:
                        line2print += "1"
                else:
                        line2print += "0"
#        print(len(line2print))
        if len(i)>10:
                print(i[0:10]+"       "+line2print)
        else:
                print(i+"       "+line2print)
