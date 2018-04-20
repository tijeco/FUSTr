import sys

from itertools import groupby

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        headerStr = header.__next__()[1:].strip()#.split("_")[0]#Entire line, add .split[0] for just first column
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)

seq_length=0
def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-i" in sys.argv:
    aln_file = getOptionValue("-i").strip("/")
else:

    print("\nplease specify input directory name using -i <directory_name> \n")
    sys.exit()
    # print(output[currentFile],input[currentFile])


with open(aln_file.split("aln")[0]+"phy","w") as out:
    sequence_iterator = fasta_iter(aln_file)
    first_line =True
    for ff in sequence_iterator:

        headerStr, seq = ff
        if first_line:
            seq_length = len(seq)
            num_lines = num_lines = sum(1 for line in open(aln_file) if line[0]=='>')
            out.write(str(num_lines)+" "+str(seq_length)+"\n")
            first_line=False

        seq_length = len(seq)
        out.write(headerStr.strip('>')+"\t")
        out.write(seq +"\n")
