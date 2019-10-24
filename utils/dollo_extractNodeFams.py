import sys
import os.path
from itertools import groupby

def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue

if "-i" in sys.argv:
    working_dir = getOptionValue("-i").strip("/")
else:
    print("\nplease specify input file name using -i <file_name> \n")
    sys.exit()

if "-o" in sys.argv:
    working_dir = getOptionValue("-o").strip("/")
else:
    print("\nplease specify output file name using -o <file_name> \n")
    sys.exit()


start = False
start_dict = {}
start_key = ""
with open(sys.argv[1]) as f:
        for line in f:
                if "Any Steps?" in line:
                        start = True

                if start:
                        #print(len(line.strip()),line.replace("                           ",''))
                        if "                           " not in line:
                                try:
                                        start_key = line.split()[0] +"_" +line.split()[1]
                                        start_dict[start_key] = line.strip().split("yes")[1]
                                except:
                                        continue
                        else:
                                try:
                                        start_dict[start_key]+=line.strip()
                                except:
                                        continue
                        #print(start_key)

#print(start_dict.keys())
for node in start_dict:
#       print(node)
#       print(len(start_dict[node].replace(" ",'')))
        for i in range(len(start_dict[node].replace(" ",''))):
                if start_dict[node].replace(" ",'')[i] == "1":
                        out.write(str(node)+" FUSTrFam_"+str(i+1))
