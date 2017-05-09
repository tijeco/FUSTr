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
    with open(working_dir+"/ChiSq.txt","w") as out:
        with open(finalStatsfile) as f:
            for line in f:
                row = line.strip().split()
                ChiSq_dict[row[0]] = True
                # try:
                #
                #     ChiSq_dict[row[0]][row[1]] = (row[2],row[3])
                # except:
                #     print(row[0])
                #     ChiSq_dict[row[0]] = {}
                #     ChiSq_dict[row[0]][row[1]] = (row[2],row[3])

print(ChiSq_dict)
