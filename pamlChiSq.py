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
                #ChiSq_dict[row[0]] = True
                try:
                    print(row[0])
                    if row[0] not in ChiSq_dict:
                        ChiSq_dict[row[0]] = {}
                    ChiSq_dict[row[0]][row[1]] = (float(row[2]),float(row[3]))
                except:
                    None
                # try:
                #
                #     ChiSq_dict[row[0]][row[1]] = (row[2],row[3])
                # except:
                #     print(row[0])
                #     ChiSq_dict[row[0]] = {}
                #     ChiSq_dict[row[0]][row[1]] = (row[2],row[3])

print(ChiSq_dict)

for i in ChiSq_dict.keys():


    print("M3-M0",i,2*(ChiSq_dict[i]["M3"][1]-ChiSq_dict[i]["M0"][1]),ChiSq_dict[i]["M3"][0]-ChiSq_dict[i]["M0"][0])
    M3.M0_chiSq = 2*(ChiSq_dict[i]["M3"][1]-ChiSq_dict[i]["M0"][1])
    M3.M0_df = ChiSq_dict[i]["M3"][0]-ChiSq_dict[i]["M0"][0]


    print("M2-M1",i,2*(ChiSq_dict[i]["M2"][1]-ChiSq_dict[i]["M1"][1]),ChiSq_dict[i]["M2"][0]-ChiSq_dict[i]["M1"][0])
    M2.M1_chiSq = 2*(ChiSq_dict[i]["M2"][1]-ChiSq_dict[i]["M1"][1])
    M2.M1_df = ChiSq_dict[i]["M2"][0]-ChiSq_dict[i]["M1"][0]


    print("M8-M7",i,2*(ChiSq_dict[i]["M8"][1]-ChiSq_dict[i]["M7"][1]),ChiSq_dict[i]["M8"][0]-ChiSq_dict[i]["M7"][0])
    M8.M7_chiSq = 2*(ChiSq_dict[i]["M8"][1]-ChiSq_dict[i]["M7"][1])
    M8.M7_df = ChiSq_dict[i]["M8"][0]-ChiSq_dict[i]["M7"][0]


    print("M8 M8a",i,2*(ChiSq_dict[i]["M8"][1]-ChiSq_dict[i]["M8a"][1]),    ChiSq_dict[i]["M8"][0]-ChiSq_dict[i]["M8a"][0])
    M8.M8a_chiSq = 2*(ChiSq_dict[i]["M8"][1]-ChiSq_dict[i]["M8a"][1])
    M8.M8a_df = ChiSq_dict[i]["M8"][0]-ChiSq_dict[i]["M8a"][0]  


    # print(i)
    # print(ChiSq_dict[i])
