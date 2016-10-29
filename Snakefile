SAMPLES = ['Sample1', 'Sample2']

rule all:
    input:
        expand('{sample}.txt', sample=SAMPLES)

rule quantify_genes:
    input:
        genome = 'genome.fa',
        r1 = 'fastq/{sample}.R1.fastq.gz',
        r2 = 'fastq/{sample}.R2.fastq.gz'
    output:
        '{sample}.txt'
    shell:
        'echo {input.genome} {input.r1} {input.r2} > {output}'
