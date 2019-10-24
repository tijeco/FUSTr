import sys
import os.path
from itertools import groupby
import pandas as pd

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-d" in sys.argv:
    working_dir = getOptionValue("-d").strip("/")
else:

    print("\nplease specify input directory name using -d <directory_name> \n")
    sys.exit()

if "-o" in sys.argv:
    out_file = getOptionValue("-o").strip("/")
else:

    print("\nplease specify output file using -o <directory_name> \n")
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
                if num not in species_dict[species]:
                    species_dict[species][num] = 0
                species_dict[species][num] += 1

pd.DataFrame.from_dict(species_dict).to_csv(out_file)
