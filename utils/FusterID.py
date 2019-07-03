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
if "-i" in sys.argv:
    input_file = getOptionValue("-i").strip("/")
else:

    print("\nplease specify input directory name using -i <directory_name> \n")
    sys.exit()



fusterID_file = working_dir+"/intermediate_files/fusterID.txt"
id_dict = {}
if os.path.exists(fusterID_file):

    #make id_dict
    with open(fusterID_file) as f:
        for line in f:
            row = line.strip().split()
            id_dict[row[0]] = row[1]


# print(testString[:len("fusterID_")])
# print(testString.split("fusterID_"))
# print([int(s) for s in testString if s.isdigit()])

def findIDs(s):
    s_list = s.split("fusterID_")
    new_string = s_list[0]
    if len(s_list) > 1:
        for i in range(len(s_list)):
            id_num = ""
            for char in s_list[i]:
                if char.isdigit():
                    id_num += char
                else:
                    break
            if id_num != "":
                if "fusterID_"+id_num in id_dict:
                    new_string += id_dict["fusterID_" + id_num] + s_list[i][len(id_num):]
                else:
                    print("ERROR: fusterID_" + id_num +" not found, exiting")
                    sys.exit()
    else:
        new_string = s
    return new_string


with open(input_file) as f:
    for line in f:
        print(findIDs(line.strip()))
