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


fusterID_file = working_dir+"/intermediate_files/fusterID.txt"
family_file = working_dir + "/final_results/famsUnderSelection.txt"
familyDir = working_dir + "/Families/"
selectDir = working_dir + "/final_results/famsUnderSelection_dir/"
os.makedirs( selectDir,  exist_ok=True)

id_dict = {}
# fams = {}
if os.path.exists(fusterID_file):
    if os.path.exists(family_file):
        #make id_dict
        with open(fusterID_file) as f:
            for line in f:
                row = line.strip().split()
                id_dict[row[0]] = row[1]
        with open(family_file) as f:
            line1 = True
            for line in f:
                if line1:
                    line1 = False
                else:
                    row = line.strip().split()
                    # fams[row[0]] = True
                    current_fam = row[0]

                    pep = familyDir + current_fam + ".fa"
                    phy = familyDir + current_fam + ".fa"
                    codon = familyDir + current_fam + "_dir/" + current_fam + ".codon.phylip"
                    tree = familyDir + current_fam + "_dir/" + current_fam + ".tree"

                    new_pep = selectDir+ pep.strip("/")[-1]
                    new_phy = selectDir+ phy.strip("/")[-1]
                    new_codon = selectDir+ codon.strip("/")[-1]
                    new_tree = selectDir+ tree.strip("/")[-1]


                    for old,new  in [(pep,new_pep),(phy,new_phy),(codon,new_codon),(tree,new_tree)]:
                        with open(new,"w") as out:
                            with open(old) as f0:
                                for line0 in f0:
                                    line2print = line0
                                    for key in id_dict:
                                        line2print = line2print.replace(key,id_dict[key])
                                    out.write(line2print)






    else:
        print(family_file+" does not exist\nexiting")
        sys.exit()
else:
    print(fusterID_file +" does not exist\nexiting")
    sys.exit()
