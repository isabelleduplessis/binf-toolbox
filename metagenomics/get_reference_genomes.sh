#!/bin/bash

# Dependencies: ncbi datasets command line tools

# Run script example: ./get_reference_genomes.sh --taxid 1313 --name s.pneumoniae --type bac

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --taxid) taxid="$2"; shift ;;
        --name) org="$2"; shift ;;
	--type) type="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p ${type}
cd ${type}

datasets download genome taxon $taxid --reference --assembly-source RefSeq --assembly-level complete --include genome,gff3 --filename ${org}.zip

unzip ${org}.zip -d ${org} || exit 1
rm ${org}.zip

mv ${org}/ncbi_dataset/data/GCF*/* ${org}/
rm -rf ${org}/ncbi_dataset
rm ${org}/md5sum.txt ${org}/README.md


f=$(basename ${org}/*.fna)
mid=$(echo "$f" | sed -E 's/^GCF_[0-9]+\.[0-9]+_([^_].*)_genomic\.fna/\1/')

mv "${org}/$f" "${org}/${org}_${mid}.fna"
mv "${org}/genomic.gff" "${org}/${org}_${mid}.gff"
