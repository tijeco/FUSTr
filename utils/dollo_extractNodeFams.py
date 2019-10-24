import sys
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
                        print(node,"FUSTrFam_"+str(i+1))
