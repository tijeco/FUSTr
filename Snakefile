rule complex_conversion:
    input:
        "{dataset}/fasta"
    output:
        "{dataset}/txt"
    shell:
        "echo {input} > {output}"
