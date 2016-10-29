from os.path import join

# Globals ---------------------------------------------------------------------

# Full path to a FASTA file.
#GENOME = 'genome.fa'

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = 'fastq/'

SAMPLES, = glob_wildcards('{sample}.fasta')
print(SAMPLES)
print(glob_wildcards('{sample}.fasta'))
print('{sample}.fasta')
#
# AMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,Samp[^/]+}.R1.fastq.gz'))
# print(AMPLES)
# print(glob_wildcards(join(FASTQ_DIR, '{sample,Samp[^/]+}.R1.fastq.gz')))

# print('{sample,Samp[^/]+}.R1.fastq.gz')

# PATTERN_R1 = '{sample}.R1.fastq.gz'
# PATTERN_R2 = '{sample}.R2.fastq.gz'

# Rules -----------------------------------------------------------------------
#
# rule all:
#     input:
#         'test.txt'
#
# rule quantify_genes:
#     input:
#         genome = '{sample}.fasta'
#         #genome = GENOME,
#         # r1 = join(FASTQ_DIR, PATTERN_R1),
#         # r2 = join(FASTQ_DIR, PATTERN_R2)
#     output:
#         '{sample}.txt'
#     shell:
#         'echo {input.genome}  > {output}'
#
# rule collate_outputs:
#     input:
#         expand('{sample}.txt', sample=SAMPLES)
#     output:
#         'test.txt'
#     run:
#         with open(output[0], 'w') as out:
#             for i in input:
#                 sample = i.split('.')[0]
#                 for line in open(i):
#                     out.write(sample + ' ' + line)
rule thingy:
    input:
        "{sample}.fasta"
    output:
        "{sample}.texas"
    shell:
        'echo {input}>output'
