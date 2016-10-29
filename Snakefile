rule complex_conversion:
    input:
        "{dataset}/fasta"
    output:
        "{dataset}/file.{group}.txt"
    shell:
        "echo --group {wildcards.group} < {input} > {output}"
