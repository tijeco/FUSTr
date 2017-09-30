import json
from pprint import pprint

with open('family_853.aln.codon.FUBAR.json') as data_file:
    data = json.load(data_file)

pprint(data["MLE"]["content"]['0'][0])
pprint(data["MLE"]["headers"])
pprint(data["MLE"]["headers"][0])
