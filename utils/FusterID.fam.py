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
    working_file = getOptionValue("-i").strip("/")
else:

    print("\nplease specify input directory name using -i <directory_name> \n")
    sys.exit()
# SAMPLES, = glob_wildcards("{sample}.fasta")

fusterID_file = working_dir+"/Temp/fusterID.txt"
family_file = working_dir + "/Temp/all.pep.combined_r90_SLX.fnodes"


id_dict = {}
if os.path.exists(fusterID_file):
    if os.path.exists(family_file):
        #make id_dict
        with open(fusterID_file) as f:
            for line in f:
                row = line.strip().split()
                id_dict[row[0]] = row[1]
        with open(working_file) as f:
            for line in f:
                row = line.strip().split()
                if row[0] in id_dict:
                    print(id_dict[row[0]]+"\t"+row[1])
                else:
                    print(line.strip())
        # with open(working_dir+"_familyFile.txt","w") as out:
        #     with open(family_file) as f:
        #         for line in f:
        #             row = line.strip().split()
        #             out.write(row[0]+"\t"+id_dict[row[1]]+"\n")
    else:
        print(family_file+" does not exist\nexiting")
        sys.exit()
else:
    print(fusterID_file +" does not exist\nexiting")
    sys.exit()
