#!/bin/bash
# nf_haplotypecaller 0.0.1

main() {

    echo "Value of sampleID: '$sampleID'"
    echo "Value of fastqR1: '$fastqR1'"
    echo "Value of fastqR2: '$fastqR2'"
    echo "Value of ref: '$ref'"
    echo "Value of fastaIn: '$fastaIn'"
    echo "Value of knownVariants: '${knownVariants[@]}'"
    echo "Value of dbsnp: '$dbsnp'"

    dx-download-all-inputs
    mkdir -p ./knownVariants
    mkdir -p ./outputs
    mkdir -p ./out/vcf
    mkdir -p ./out/tbi
    mkdir -p ./out/intermediate_tar
    mkdir -p out/report
    mkdir -p out/dag
    for i in ${knownVariants_path[@]}
    do
        mv ${i} ./knownVariants
    done

    nextflow run main.nf \
        --sampleID ${sampleID} \
        --fastqR1 ${fastqR1_path} \
        --fastqR2 ${fastqR2_path} \
        --ref ${ref_path} \
        --fastaIn ${fastaIn_path} \
        --knownVariants ./knownVariants/* \
        --dbsnp ${dbsnp_path} \
        --outDir ./outputs \
        -with-trace false -with-report report.html -with-timeline false -with-dag dag.png

    cp ~/outputs/${sampleID}.vcf.gz ./out/vcf/
    cp ~/outputs/${sampleID}.vcf.gz.tbi ./out/tbi/
    tar czf ~/out/intermediate_tar/output.tar.gz ~/outputs
    mv *.html ~/out/report/
    mv *.png ~/out/dag/
    dx-upload-all-outputs
}
