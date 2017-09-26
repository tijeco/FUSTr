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
# if "-i" in sys.argv:
#     working_file = getOptionValue("-i").strip("/")
# else:
#
#     print("\nplease specify input directory name using -i <directory_name> \n")
#     sys.exit()
# SAMPLES, = glob_wildcards("{sample}.fasta")

fusterID_file = working_dir+"/intermediate_files/fusterID.txt"
blast_file = working_dir + "/intermediate_files/all.pep.combined.blastall.out"


id_dict = {}
if os.path.exists(fusterID_file):
    if os.path.exists(blast_file):
        #make id_dict
        with open(fusterID_file) as f:
            for line in f:
                row = line.strip().split()
                id_dict[row[0]] = row[1]
        with open(blast_file) as f:
            for line in f:
                row = line.strip().split()
                line2print = ""
                for i in row:
                    if i in id_dict:
                        line2print+=id_dict[i]+" "
                    else:
                        line2print+= i +" "
                print(line2print)
