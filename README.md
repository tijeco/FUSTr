# FUSTr
Families Under Selection in Transcriptomes
# Introduction  
Fuster is a pipeline that clusters coding sequences from transcriptomes into protein families, and then analyzes those families for positive selection.

# Getting started

The only software needed to run FUSTr is [Docker](https://www.docker.com/). FUSTr takes as input a directory containing transcriptome assemblies

**Important:** transcriptome assemblies need to be fasta files ending in .fasta, all in one directory.


Download FUSTr with the following command
```bash
git clone https://github.com/tijeco/FUSTr.git
```


With Docker installed correctly on your system issue the following command to initialize the Docker container

```
bash FUSTr/setup_docker.sh <directory_containing_fastas>
```

Once the docker container has been initialized you can enter it using the following command

```
docker run -it fustr /bin/bash
```

Now that you are in the docker container, your data is in /home/usr/data, to run FUSTr simply issue the following command
```
FUSTr -d ./data -t <number_of_threads>
```

The output will be in data/file

You can use scp to transfer this to your local machine.
