configfile: "config/config.yaml"

import glob

# list/comment-out analyses output (needs to folder for some, file for others)

analyses=[
    "flye",
    "medaka",
    "checkm2",
    "decon",
    "restart",
    "checkm2_decon",
    "gtdbtk",
    "genomad",
    "bakta",
]

#-----------------------------------------------------------------------------#
# some helper functions to gather variable number of read files and sample ids
def glob_ids():
    fqs=glob.glob("analysis/reads/*.fq.gz") 
    return sorted({
        re.sub("-.*.fq.gz$", '', fq.removeprefix("analysis/reads/")) for fq in fqs
    })

def glob_fqs(wildcards):
    return glob.glob(f"analysis/reads/{wildcards.x}-*.fq.gz")


#-----------------------------------------------------------------------------#
FQP, = glob_wildcards("analysis/reads/{z}.fq.gz")
rule all:
    input:
        expand("analysis/sylph/{z}.tsv", z=FQP),
        expand("analysis/{a}/{x}", a=analyses, x=glob_ids())

# Taxonomic profiling of reads to check for off samples and contaminations
# run sample by sample so it's easier to recompute for added samples
rule sylph:
    conda: "envs/sylph.yaml"
    threads: 4
    input: "analysis/reads/{z}.fq.gz"
    output:
        tsv="analysis/sylph/{z}.tsv"
    shell:
        "sylph profile -t {threads} {config[sylph_db]} {input} > {output.tsv}"

# Assembly of reads and postprocessing for consistent contig naming
rule flye:
    conda: "envs/flye.yaml"
    input: glob_fqs
    output:
        dir=directory("analysis/flye/{x}/"),
        raw="analysis/flye/{x}/assembly.fasta",
        tsv="analysis/flye/{x}/{x}.tsv",
        fna="analysis/flye/{x}/{x}.fna"
    threads: workflow.cores
    shell:
        "flye -t {threads} {config[flye]} --nano-hq {input} -o {output.dir}; "
        "flye-post {config[flye_post]} {output.dir}"

# medaka needs a single fq library > temporary merge fqs
# medaka can screw with sorting order of contigs > reorder
rule medaka:
    conda: "envs/medaka.yaml"
    threads: 4
    input:
        fqz=glob_fqs,
        fna="analysis/flye/{x}/{x}.fna"
    output:
        dir=directory("analysis/medaka/{x}"),
        tmp=temp("analysis/medaka/{x}/tmp.fq"),
        cns="analysis/medaka/{x}/consensus.fasta",
        fna="analysis/medaka/{x}/{x}.fna"
    shell:
        "seqkit seq -o {output.tmp} {input.fqz}; "
        "medaka_consensus -t {threads} -o {output.dir} -i {output.tmp}"
        " -d {input.fna};"
        "seqkit sort -rl {output.cns} | seqkit seq -o {output.fna};"

# check completeness and potential contamination using marker genes
rule checkm2:
    conda: "envs/checkm2.yaml"
    threads: 4
    input:
        "analysis/medaka/{x}/{x}.fna"
    output:
        dir=directory("analysis/checkm2/{x}/"),
        tsv="analysis/checkm2/{x}/quality_report.tsv"
    shell:
        "checkm2 predict -t {threads} --force"
        " --database_path {config[checkm2_db]} -i {input} -o {output.dir};"

rule checkm2_collect:
    conda: "envs/decon.yaml"
    threads: 1
    input: expand("analysis/checkm2/{x}/quality_report.tsv", x=glob_ids())
    output:
        all="analysis/checkm2/checkm2_all.tsv",
        bad="analysis/checkm2/checkm2_contaminated.tsv"
    shell:    
        "csvtk -tT concat {input}.tsv > {output.all}; "
        "csvtk -tT filter -f 'Contamination>5' {output.all} > {output.bad}"

# remove low-coverage contaminations from assemblies flagged by checkm2
rule decon:
    conda: "envs/decon.yaml"
    threads: 1
    input:
        fna="analysis/medaka/{x}/{x}.fna",
        tsv="analysis/flye/{x}/{x}.tsv",
        ckm="analysis/checkm2/{x}/quality_report.tsv"
    output:
        dir=directory("analysis/decon/{x}/"),
        fna="analysis/decon/{x}/{x}.fna",
        tsv="analysis/decon/{x}/{x}.tsv"
    shell:
        """
        set -euo pipefail
        # checks checkm2 third field ("Contamination") of second row (ignoring header)
        is_contaminated() {{ # <checkm2.tsv> <threshold>
            # NOTE: force numerical comparison on awk with "X"+0
            awk -F'\\t' -v t="$2" 'NR==2 {{exit ($3+0 > t+0) ? 0 : 1}}' "$1"
        }}

        if is_contaminated {input.ckm} {config[decon][min_contamination]}; then
        # filter assembly fasta and tsv
          csvtk -tT filter -f 'coverage>={config[decon][min_coverage]}' {input.tsv} > {output.tsv}
          seqkit grep -rf <(cut -f1 {output.tsv}) {input.fna} > {output.fna}
        else
          # nothing to do - symlink
          scripts/symlink.sh {input.fna} {output.fna}
          scripts/symlink.sh {input.tsv} {output.tsv}
        fi
        """

# start/orient circular main chromosomes on dnaN gene
rule restart:
    input:
        fna="analysis/decon/{x}/{x}.fna",
        tsv="analysis/decon/{x}/{x}.tsv"
    output:
        dir=directory("analysis/restart/{x}/"),
        fna="analysis/restart/{x}/{x}.fna",
        tsv="analysis/restart/{x}/{x}.tsv"
    shell:
        """
        # checks flye_info.tsv 4th field ("circular") of second row/first contig (ignoring header)
        is_circular() {{
            awk -F'\\t' 'NR==2 {{exit ($4 == "Y") ? 0 : 1}}' "$1"  
        }}

        if is_circular {input.tsv}; then
          scripts/seq-circ-restart --status closed --start {config[seq-circ-restart_db]} -o {output.fna} {input.fna}
          scripts/symlink.sh {input.tsv} {output.tsv}
        else
          scripts/symlink.sh {input.fna} {output.fna}
          scripts/symlink.sh {input.tsv} {output.tsv}
        fi
        """

rule checkm2_decon:
    conda: "envs/checkm2.yaml"
    threads: 4
    input:
        "analysis/restart/{x}/{x}.fna"
    output:
        dir=directory("analysis/checkm2_decon/{x}/"),
        tsv="analysis/checkm2_decon/{x}/quality_report.tsv"
    shell:
        "checkm2 predict -t {threads} --force"
        " --database_path {config[checkm2_db]} -i {input} -o {output.dir};"

rule checkm2_decon_collect:
    conda: "envs/decon.yaml"
    threads: 1
    input: expand("analysis/checkm2_decon/{x}/quality_report.tsv", x=glob_ids())
    output:
        all="analysis/checkm2_decon/checkm2_decon_all.tsv",
        bad="analysis/checkm2_decon/checkm2_decon_contaminated.tsv"
    shell:
        """
        csvtk -tT concat {input} > {output.all}
        csvtk -tT filter -f 'Contamination>5' {output.all} > {output.bad}
        """

# taxonomically classify genomes using GTDB
rule gtdbtk:
    conda: "envs/gtdbtk.yaml"
    threads: 4
    input:
        dir="analysis/restart/{x}/"
    output:
        dir=directory("analysis/gtdbtk/{x}/"),
        gtd="analysis/gtdbtk/{x}/gtdbtk.bac120.summary.tsv"
    shell:
        "export GTDBTK_DATA_PATH={config[gtdbtk_db]}; "
        "gtdbtk classify_wf --cpus {threads} --genome_dir {input.dir}"
        " --out_dir {output.dir} --mash_db {config[gtdbtk_db]}"

rule gtdbtk_collect:
    conda: "envs/decon.yaml"
    threads: 1
    input: expand("analysis/gtdbtk/{x}/gtdbtk.bac120.summary.tsv", x=glob_ids())
    output:
        all="analysis/gtdbtk/gtdbtk_bac120_all.tsv"
    shell:    
        "csvtk -tT concat {input} > {output.all}"

# identify plasmids and prophages
rule genomad:
    conda: "envs/genomad.yaml"
    threads: 4
    input:
        fna="analysis/restart/{x}/{x}.fna"
    output:
        dir=directory("analysis/genomad/{x}"),
        pls="analysis/genomad/{x}/{x}_summary/{x}_plasmid_summary.tsv"
    shell:
        "genomad end-to-end --cleanup {input.fna} {output.dir} {config[genomad_db]}"

# annotate genomes
rule bakta:
    conda: "envs/bakta.yaml"
    threads: workflow.cores
    input:
        fna="analysis/restart/{x}/{x}.fna",
        tsv="analysis/restart/{x}/{x}.tsv",
        pls="analysis/genomad/{x}/{x}_summary/{x}_plasmid_summary.tsv",
        gtd="analysis/gtdbtk/{x}/gtdbtk.bac120.summary.tsv"
    output:
        dir=directory("analysis/bakta/{x}/"),
        rpl="analysis/bakta/{x}-repl.tsv"
    shell:
        """
        Rscript --vanilla scripts/bakta-replicons.R {input.tsv} {input.pls} {output.rpl}

        x={wildcards.x}
        gram="+"
        if [ $(grep -w vup2 "databases/gram-negative.tsv") ]; then gram="-"; fi;
        
        species=$(tail -1 {input.gtd} | cut -f2 | sed 's/.*;s__//')
        # capture unclassified genus_species in GTDB
        if [ -z "$species"]; then
          species="Incertae sedis"
        fi;
        
        genus=${{species%% *}}
        species=${{species#* }}

        bakta -t {threads} --db {config[bakta_db]} --prefix $x \
          --output {output.dir} --genus $genus --species $species --strain $x \
          --gram $gram --locus-tag ${{x^^}} --keep-contig-headers --compliant \
         --replicons {output.rpl} {input.fna}
        """

    
