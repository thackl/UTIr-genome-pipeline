configfile: "config/config.yaml"

import glob

# list/comment-out analyses results (needs to folder for some, file for others)
analyses=[
    "sylph",
    "flye",
    "checkm2",
    "gtdbtk",
    # "genomad",
    # "bakta",
]

#-----------------------------------------------------------------------------#
# some helper functions to gather variable number of read files and sample ids
def glob_ids():
    fqs=glob.glob("input/*.fq.gz") 
    return sorted({
        re.sub("-.*.fq.gz$", '', fq.removeprefix("input/")) for fq in fqs
    })

def glob_fqs(wildcards):
    return glob.glob(f"input/{wildcards.x}-*.fq.gz")

print(glob_ids())

#-----------------------------------------------------------------------------#
rule all:
    input:
        expand("results/{x}/{y}", x=glob_ids(), y=analyses)

rule assemble_flye:
    conda: "envs/flye.yaml"
    input: glob_fqs
    output:
        fna="results/{x}/flye/assembly.fasta",
        dir=directory("results/{x}/flye")
    threads: workflow.cores
    shell:
        "flye -t {threads} {config[flye]} --nano-hq {input} -o {output.dir};"

rule assemble_flye_post:
    input: "results/{x}/flye/assembly.fasta"
    output:
        fna="results/{x}/flye/{x}.fna",
        tsv="results/{x}/flye/{x}.tsv"
    threads: 1
    shell:
        "flye-post {config[flye_post]} {input} {wildcards.x}"

# medaka needs a single fq library > temporary merge fqs
# medaka can screw with sorting order of contigs > reorder
rule polish_medaka:
    conda: "envs/medaka.yaml"
    threads: 4
    input:
        fqz=glob_fqs,
        fna="results/{x}/flye/{x}.fna"
    output:
        dir=directory("results/{x}/medaka/"),
        tmp=temp("results/{x}/medaka/tmp.fq"),
        cns="results/{x}/medaka/consensus.fasta",
        fna="results/{x}/medaka/{x}.fna"
    shell:
        "seqkit seq -o {output.tmp} {input.fqz}; "
        "medaka_consensus -t {threads} -o {output.dir} -i {output.tmp}"
        " -d {input.fna};"
        "seqkit sort -rl {output.cns} | seqkit seq -o {output.fna};"

rule assess_checkm2:
    conda: "envs/checkm2.yaml"
    threads: 4
    input:
        "results/{x}/medaka/{x}.fna"
        #"results/{x}/flye/{x}.fna"
    output:
        dir=directory("results/{x}/checkm2/"),
        tsv="results/{x}/checkm2/quality_report.tsv"
    shell:
        "checkm2 predict -t {threads} --force"
        " --database_path {config[checkm2_db]} -i {input} -o {output.dir};"

rule taxify_sylph:
    conda: "envs/sylph.yaml"
    threads: 4
    input: glob_fqs
    output:
        dir=directory("results/{x}/sylph"),
        tsv="results/{x}/sylph/{x}-sylph-profiling.tsv"
    shell:
        "sylph profile -t {threads} {config[sylph_db]} {input} > {output.tsv}"

rule taxify_gtdbtk:
    conda: "envs/gtdbtk.yaml"
    threads: 4
    input:
        dir="results/{x}/medaka/"
    output:
        dir=directory("results/{x}/gtdbtk")
    shell:
        "export GTDBTK_DATA_PATH={config[gtdbtk_db]}; "
        "gtdbtk classify_wf --cpus {threads} --genome_dir {input.dir}"
        " --out_dir {output.dir} --mash_db {config[gtdbtk_db]}"
        
rule annotate_bakta:
    conda: "envs/bakta.yaml"
    threads: workflow.cores
    input:
        fna="results/{x}/medaka/{x}.fna"
    output:
        dir=directory("results/{x}/bakta/")
    shell:
        "bakta -t {threads} --db {config[bakta_db]} --output {output.dir}"
        " --prefix {wildcards.x} --locus {wildcards.x} {input.fna}"

rule annotate_genomad:
    conda: "envs/genomad.yaml"
    threads: 4
    input:
        fna="results/{x}/medaka/{x}.fna"
    output:
        dir=directory("results/{x}/genomad/")
    shell:
        "genomad end-to-end --cleanup {input.fna} {output.dir} {config[genomad_db]}"
    
