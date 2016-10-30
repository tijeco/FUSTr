#!/usr/bin/python3
from itertools import groupby


def fasta_iter(fasta_name):


    fh = open(fasta_name)


    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        headerStr = header.__next__()[1:].strip()
        # print(header)


        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

fiter = fasta_iter('1.fasta')

for ff in fiter:

    headerStr, seq = ff

    print(">"+headerStr)
    print(seq)
    print(str(3))
    print("a" in "abc")
