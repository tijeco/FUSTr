import sys
import os.path


def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-d" in sys.argv:
    working_dir = getOptionValue("-d").strip("/")
else:

    print("\nplease specify input directory name using -d <directory_name> \n")
    sys.exit()

finalStatsfile = working_dir + "/finalStatsfile.txt"

ChiSq_dict = {}
if os.path.exists(finalStatsfile):
    with open(working_dir.split("/")[1]+"_ChiSq.txt") as out:
        with open(finalStatsfile) as f:
            for line in f:
                row = line.strip().split()
                if row[0] not in ChiSq_dict:
                    ChiSq_dict[row[0]] = {}
                ChiSq_dict[row[0]][row[1]] = (row[2],row[3])


print(ChiSq_dict)
