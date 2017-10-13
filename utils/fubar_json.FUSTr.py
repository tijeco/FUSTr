import os
import json
from pprint import pprint
import sys

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
def srcdir(str):
    fullPath = os.path.dirname(os.path.abspath(str.split("/")[-1]))# +"/" +str
    return fullPath
doAll=False
if "-d" in sys.argv:
    input_directory = getOptionValue("-d")

elif "-all" in sys.argv:
    json_file = getOptionValue("-all")
    doAll=True
else:

    # input_directory = "data/"
    print("\nplease specify input directory name using -d <file_name> \n")
    sys.exit()

if doAll:
    csv_file = json_file.strip("json")+"csvv"
    os.makedirs("all_csv", exist_ok=True)
    with open(csv_file,"w") as out:
        with open(json_file) as data_file:
            data = json.load(data_file)

        line2print=""
        for i in range(len(data["MLE"]["headers"])):
            line2print+=data["MLE"]["headers"][i][0]+" " +data["MLE"]["headers"][i][1]+","
        out.write(line2print[:-1]+'\n')

        for row in data["MLE"]["content"]['0']:
            line2print=""
            for i in row:
                line2print += str(i)+","
            out.write(line2print[:-1]+'\n')




else:

    fubar_file = input_directory.strip('/')  +"/final_results/famsUnderSelection.txt"

    famsUnderSelectionDict = {}
    line1 = True
    with open(fubar_file) as f:
        for line in f:
            if line1:
                line1 = False
                continue
            row = line.strip().split()
            famsUnderSelectionDict[row[0]] = True

    # print(famsUnderSelectionDict)
    os.makedirs(input_directory.strip('/')+"/fubar_out", exist_ok=True)
    for i in famsUnderSelectionDict.keys():
        # print(input_directory.strip('/')+"/Families/"+i+"_dir/"+i+".aln.codon.FUBAR.json")
        current_json = input_directory.strip('/')+"/Families/"+i+"_dir/"+i+".aln.codon.FUBAR.json"

        current_csv = input_directory.strip('/')+"/fubar_out/"+i+".fubar.csv"
        print("writing to "+current_csv)
        with open(current_csv,"w") as out:
            with open(current_json) as data_file:
                data = json.load(data_file)

            line2print=""
            for i in range(len(data["MLE"]["headers"])):
                line2print+=data["MLE"]["headers"][i][0]+" " +data["MLE"]["headers"][i][1]+","
            out.write(line2print[:-1]+'\n')

            for row in data["MLE"]["content"]['0']:
                line2print=""
                for i in row:
                    line2print += str(i)+","
                out.write(line2print[:-1]+'\n')
