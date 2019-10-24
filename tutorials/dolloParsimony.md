# Gene family loss/gain analysis with dollo parsimony using families from FUSTr

## Requirements

1. ```FUSTr.fams``` file
1. Species delimiting character
  2. within the gene name I am assuming I can get a unique identifier for the species by extracting the first field of each gene name  using "_"
  3. dollop requires species names only be 10 characters, so typically this works best if only the genus is kept
  4. Your data may need to be reconfigured if it doesn't meet this assumption
1. Species tree
1. dollop software





## convert ```FUSTr.fams``` file to binary matrix for dollop


```bash
python FUSTr/utils/FUSTrFams2DolloMatrix.py -d path/to/FUSTr_input_directory >  ~/infile
```

I have suggested the output be written to ```~/infile```, this is because dollop works best when the input file is named this way and can be found in the home directory

## species tree
For dollop, the species tree needs to be saved as ```~/intree```

## dollop

dollop is one of those command line utilities where you are prompted with a screen to choose some options. It is best to launch dollop from your home directory

```bash
dollop
```
The options you need to change are the following
* Search for best tree?
  * Type "U" to change this
* Print out steps in each character
  * Type "4" to change this to "Yes"
* Print states at all nodes of tree
  * Type "5" to change this to "Yes"

This should generate two files, ```outfile``` and ```outtree```

The outtree will have node labels that will be used in the next two steps
# determine how many gene families were lost/gained at each node
```bash
python parseDollo.py outfile > node.loss_gain.csv
```
# determine which specific families were gained at each node

```bash
python dollo_extractNodeFams.py -i outfile -o node_fams.txt
```
