import sys
from itertools import groupby


def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # Entire line, add .split[0] for just first column
        headerStr = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)


untrimmed = sys.argv[1]
column_file = sys.argv[2]
nucleotide = sys.argv[3]

output1 = sys.argv[1].split(".")[0] + ".codon.aln"

cut = ""
longIsoform_CDS_combined = {}
sequence_iterator = fasta_iter(nucleotide)
with open(column_file) as f:
    for line in f:
        if "#ColumnsMap" in line:
            cut += line.strip().split("#ColumnsMap")[1]
    print(cut)
    cut = cut.split(',')
    cut = list(map(int, cut))
for ff in sequence_iterator:
        headerStr, seq = ff
        GeneID = headerStr
        if GeneID not in longIsoform_CDS_combined:
            longIsoform_CDS_combined[GeneID] = seq
    # Open outout
    # print(len(longIsoform_CDS_combined))
        with open(output1, "w") as out:
            # Get  column cut file

            # Get corresponding untrimmed Alignments, as original, line by line
            line1 = True
            first_line = True
            with open(untrimmed) as f:
                for line in f:
                    if line1:
                        line1 = False
                        continue

                    row = line.strip().split()
                    original = row[1]  # cds
                    header = row[0]

                    # NOTE, potetintal bug below, if exception then sequence isn't declared and it can't go forward, use continue probably
                    try:
                        sequence = longIsoform_CDS_combined[header]  # original
                    except:
                        continue
                    CodonPos = {}
                    position = 0
                    codon = ""
                    number = 1
                    for i in sequence:

                        codon += i
                        if position % 3 == 2:
                            CodonPos[number] = codon
                            number += 1
                        position += 1

                        if position % 3 == 0:
                            codon = ""
                    aaPos = 0
                    firstAA = True
                    alnPos = 0
                    prot = ""
                    trimmed = ""
                    for i in original:
                        if i != "-":
                            aaPos += 1

                        if alnPos in cut:
                            prot += i
                            if i != "-":
                                # print(aaPos,CodonPos[aaPos])
                                trimmed += CodonPos[aaPos]
                            else:
                                trimmed += "---"
                        alnPos += 1
                    num_lines = sum(1 for line in open(untrimmed))

                    out.write(">" + header + "\n")
                    out.write(trimmed + "\n")
