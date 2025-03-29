#!env bash

min_contamination=3

usage(){
cat <<EOF
Usage:
  decontaminate.sh FQ-SYLPH.TSV

  Takes a slyph profiles of individual fastq file. If it matches multiple
  genomes, maps all reads against all matched genomes and discards all reads
  matching to any genome but the primary one. Returns the cleaned fq file.

  -o  output folder
  -g  prefix to GTDB database skani folder
  -m  Minimum contamination to trigger cleaning [$min_contamination]


EOF

exit 0;
}



ARGS=`getopt --name "decontaminate.sh" \
    --options "m:o:g:h" \
    -- "$@"`

#Bad arguments
[ $? -ne 0 ] && exit 1;

# A little magic
eval set -- "$ARGS"

# Now go through all the options
while true; do
    case "$1" in
        -m) min_contamination=$2; shift 2;;
        -g) gtdb_prefix=$2; shift 2;;
        -o) outdir=$2; shift 2;;
        -h) usage && exit 0;;
        --) shift; break;;
        *) echo "$1: Unknown option" 1>&2 && exit 1;;
    esac
done

mkdir -p $outdir

for sylph_tsv in $@; do
	
    samples=( $(cut -f1 $sylph_tsv) )
    fq_in=${samples[1]}
    genomes=( $(cut -f2 $sylph_tsv | cut -f2- -d/ ) ) # loose first directory
    seq_perc=( $(cut -f4 $sylph_tsv | cut -f1 -d.) ) # need integers for bash

    genomes=("${genomes[@]/#/$gtdb_prefix}")
    n_genomes=${#genomes[@]}


    # if there are more than 2 line (header + primary genome)
    if [[ $n_genomes -gt 2 && ${seq_perc[2]} -ge 3 ]] ; then
	echo contaminated ${samples[2]};
	echo ${seq_perc[2]};
	echo ${genomes[@]};
    fi;

    bad_genomes=("${genomes[@]:2:$n_genomes}")

    genomes_fna=$outdir/genomes.fna
    genomes_paf=$outdir/genomes.paf
    fq_out=$outdir/$(basename $fq_in)
    seqkit seq ${genomes[1]} > $genomes_fna
    seqkit replace -p "^" -r DISCARD_ ${bad_genomes[@]} >> $genomes_fna

    minimap2 -c -x map-ont --secondary=no $genomes_fna $fq_in > $genomes_paf

    grep -P "\tDISCARD_" $genomes_paf | cut -f1 | seqkit grep -vf - -o $fq_out $fq_in 

    seqkit stat $fq_in $fq_out

done;
