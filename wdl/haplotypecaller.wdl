version 1.0

workflow haplotypecaller_wdl {
    input {
        String sampleID
        File fastqR1
        File fastqR2
        File ref
        File fastaIn
        Array[File] knownVariants
        File dbsnp
    }
    call mapping {
        input:
            sampleID = sampleID,
            fastqR1 = fastqR1,
            fastqR2 = fastqR2,
            ref = ref
    }
    call markdup {
        input:
            sampleID = sampleID,
            bam = mapping.bam
    }
    call fastaIndex {
        input:
            fastaIn = fastaIn
    }
    call recal {
        input:
            sampleID = sampleID,
            bam = markdup.bam_md,
            bai = markdup.bai_md,
            fasta = fastaIndex.fasta,
            fai = fastaIndex.fai,
            knownVariants = knownVariants
    }
    call bamIndex {
        input:
            bam = recal.bam_recal
    }
    call haplotypecaller {
        input:
            sampleID = sampleID,
            bam = recal.bam_recal,
            bai = bamIndex.bai_recal,
            fasta = fastaIndex.fasta,
            fai = fastaIndex.fai,
            dict = recal.dict,
            dbsnp = dbsnp
    }
    call genotype {
        input:
            sampleID = sampleID,
            fasta = fastaIndex.fasta,
            fai = fastaIndex.fai,
            dict = recal.dict,
            dbsnp = dbsnp,
            gvcf = haplotypecaller.gvcf
    }
    output {
        File vcf = genotype.vcf
        File tbi = genotype.tbi
    }
    parameter_meta {
        sampleID: {
            label: "Sample ID",
            default: "HCC1187"
        }
        fastqR1: {
            label: "Read 1 fastq (gzipped)",
            patterns: ["*.fastq.gz","*.fq.gz"],
            stream: true,
            default: "dx://wdl.vs.nf:/inputs/fastqs/subsample/HCC-1187_R1.subsampled.fastq.gz"
        }
        fastqR2: {
            label: "Read 2 fastq (gzipped)",
            patterns: ["*.fastq.gz","*.fq.gz"],
            stream: true,
            default: "dx://wdl.vs.nf:/inputs/fastqs/subsample/HCC-1187_R2.subsampled.fastq.gz"
        }
        ref: {
            label: "Reference for bwa-mem (gzipped tarball)",
            patterns: ["*.bwa-index.tar.gz"],
            suggestions: [
                {
                    project: "project-BQpp3Y804Y0xbyG4GJPQ01xv",
                    path: "/",
                    region: "aws:us-east-1",
                    name: "Reference Genome Files: AWS US (East)"
                }
            ],
            stream: true,
            default: "dx://wdl.vs.nf:/inputs/ref/hs38DH.bwa-index.tar.gz"

        }
        fastaIn: {
            label: "Reference fasta",
            patterns: ["*.fa.gz", "*.fasta.gz"],
            suggestions: [
                {
                    name: "DNAnexus Reference Genomes: AWS US-east",
                    project: "project-BQpp3Y804Y0xbyG4GJPQ01xv",
                    path: "/"
                }
            ],
            default: "dx://wdl.vs.nf:/inputs/ref/hs38DH.fa.gz"
        }
        knownVariants: {
            label: "Known variants (gzipped)",
            patterns: ["*.vcf.gz"],
            suggestions: [
                {
                    project: "project-B6JG85Z2J35vb6Z7pQ9Q02j8",
                    path: "/gatk.resources.b37",
                    name: "GRCh37/hs37d5/b37 (1000G Phase I/II) [tick one dbSNP and both indel files]: AWS US-east"
                },
                {
                    project: "project-B6JG85Z2J35vb6Z7pQ9Q02j8",
                    path: "/gatk.resources.hg19",
                    name: "UCSC hg19 [tick one dbSNP and both indel files]: AWS US-east"
                },
                {
                    project: "project-B6JG85Z2J35vb6Z7pQ9Q02j8",
                    path: "/gatk.resources.GRCh38",
                    name: "GRCh38/hg38 [tick one dbSNP and both indel files]: AWS US-east"
                }
            ]
        }
        dbsnp: {
            label: "dbSNP (gzipped)",
            patterns: ["*.vcf.gz"],
            suggestions: [
                {
                    project: "project-B6JG85Z2J35vb6Z7pQ9Q02j8",
                    path: "/gatk.resources.b37",
                    name: "GRCh37/hs37d5/b37 (1000G Phase I/II) [tick one dbSNP and both indel files]: AWS US-east"
                },
                {
                    project: "project-B6JG85Z2J35vb6Z7pQ9Q02j8",
                    path: "/gatk.resources.hg19",
                    name: "UCSC hg19 [tick one dbSNP and both indel files]: AWS US-east"
                },
                {
                    project: "project-B6JG85Z2J35vb6Z7pQ9Q02j8",
                    path: "/gatk.resources.GRCh38",
                    name: "GRCh38/hg38 [tick one dbSNP and both indel files]: AWS US-east"
                }
            ],
            default: "dx://wdl.vs.nf:/inputs/ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
        }
    }
}

task mapping {
    input {
        String sampleID
        File fastqR1
        File fastqR2
        File ref
    }
    command <<<
        mkdir -p genome
        zcat "~{ref}" | tar xvf - -C genome
        genome_file=`ls genome/*.bwt`
        genome_file="${genome_file%.bwt}"
        bwa mem -K 100000000 -t 16 -M \
            ${genome_file} \
            ~{fastqR1} ~{fastqR2} \
            -R "@RG\tID:None\tPL:None\tPU:None\tLB:None\tSM:~{sampleID}" | \
            samtools sort --threads 16 -m 2G - > ~{sampleID}.bam
        samtools index ~{sampleID}.bam
    >>>
    runtime {
        dx_instance_type: "mem2_ssd1_v2_x32"
        docker: "quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1"
    }
    output {
        File bam = "~{sampleID}.bam"
        File bai = "~{sampleID}.bam.bai"
    }
}

task markdup {
    input {
        String sampleID
        File bam
    }
    command <<<
        gatk --java-options "-Xmx40g -Xms6000m" \
            MarkDuplicates \
                --INPUT ~{bam} \
                --METRICS_FILE ~{sampleID}.md.bam.metrics \
                --TMP_DIR . \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_INDEX true \
                --OUTPUT ~{sampleID}.md.bam
        mv ~{sampleID}.md.bai ~{sampleID}.md.bam.bai
    >>>
    runtime {
        dx_instance_type: "mem3_ssd1_v2_x8"
        docker: "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }
    output {
        File bam_md = "~{sampleID}.md.bam"
        File bai_md = "~{sampleID}.md.bam.bai"
    }
}

task fastaIndex {
    input {
        File fastaIn
    }
    command <<<
        mkdir -p fasta
        cp ~{fastaIn} ./fasta/
        fasta=`ls ./fasta/*.fa.gz`
        zcat ${fasta} > ${fasta%.gz}
        fasta="${fasta%.gz}"
        samtools faidx ${fasta}
    >>>
    runtime {
        dx_instance_type: "mem1_ssd1_v2_x2"
        docker: "quay.io/biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:c17ce694dd57ab0ac1a2b86bb214e65fedef760e-0"
    }
    output {
        File fasta = "fasta/" + basename(fastaIn,".gz")
        File fai = "fasta/" + basename(fastaIn,".gz") + ".fai"
    }
}

task recal {
    input {
        String sampleID
        File bam
        File bai
        File fasta
        File fai
        Array[File] knownVariants
    }
    command <<<
        mkdir -p fasta
        cp ~{fasta} ./fasta/
        fasta=`ls ./fasta/*.fa`
        cp ~{fai} ./fasta/
        fai=`ls ./fasta/*.fai`
        gatk --java-options -Xmx40g \
            CreateSequenceDictionary \
                -R "${fasta}" \
                -O "${fasta%.fa}.dict"
        known=""
        for knownVariant in `echo ~{sep=',' knownVariants} | tr ',' '\n'`; do
            knownVariantBase=`basename "${knownVariant}"`
            zcat "${knownVariant}" > ./${knownVariantBase%.gz}
            gatk --java-options -Xmx40g \
                IndexFeatureFile \
                -I ${knownVariantBase%.gz}
            known="${known} --known-sites ./${knownVariantBase%.gz}"
        done
        gatk --java-options -Xmx40g \
            BaseRecalibrator \
                ${known} \
                -I ~{bam} \
                -O ~{sampleID}.recal.table \
                --tmp-dir . \
                -R ${fasta} \
                --use-original-qualities \
                --verbosity INFO
        gatk --java-options -Xmx40g \
            ApplyBQSR \
                -R ${fasta} \
                --input ~{bam} \
                --output ~{sampleID}.recal.bam \
                --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --use-original-qualities \
                --bqsr-recal-file ~{sampleID}.recal.table
        ls
    >>>
    runtime {
        dx_instance_type: "mem3_ssd1_v2_x8"
        docker: "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }
    output {
        File bam_recal = "~{sampleID}.recal.bam"
        File dict = "fasta/" + basename(basename(fasta,".gz"),".fa") + ".dict"
    }
}

task bamIndex {
    input {
        File bam
    }
    command <<<
        mkdir -p bam
        cp ~{bam} ./bam/
        bam=`ls ./bam/*.bam`
        samtools index ${bam}
    >>>
    runtime {
        dx_instance_type: "mem1_ssd1_x2"
        docker: "quay.io/biocontainers/mulled-v2-0560a8046fc82aa4338588eca29ff18edab2c5aa:c17ce694dd57ab0ac1a2b86bb214e65fedef760e-0"
    }
    output {
        File bai_recal = "bam/" + basename(bam) + ".bai"
    }
}

task haplotypecaller {
    input {
        String sampleID
        File bam
        File bai
        File fasta
        File fai
        File dict
        File dbsnp
    }
    command <<<
        mkdir -p fasta
        cp ~{fasta} ./fasta/
        fasta=`ls ./fasta/*.fa`
        cp ~{fai} ./fasta/
        fai=`ls ./fasta/.fai`
        cp ~{dict} ./fasta/
        dict=`ls ./fasta/*.dict`
        mkdir -p dbsnp
        cp ~{dbsnp} ./dbsnp/
        dbsnp=`ls ./dbsnp/*.vcf.gz`
        zcat "${dbsnp}" > ./${dbsnp%.gz}
        dbsnp=${dbsnp%.gz}
        gatk --java-options -Xmx40g \
            IndexFeatureFile \
            -I ${dbsnp}
        gatk --java-options "-Xmx40g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
                -R ${fasta} \
                -I ~{bam} \
                --dbsnp ${dbsnp} \
                -O ~{sampleID}.g.vcf \
                -ERC GVCF
    >>>
    runtime {
        dx_instance_type: "mem3_ssd1_v2_x8"
        docker: "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }
    output {
        File gvcf = "~{sampleID}.g.vcf"
    }
}

task genotype {
    input {
        String sampleID
        File fasta
        File fai
        File dict
        File dbsnp
        File gvcf
    }
    command <<<
        mkdir -p fasta
        cp ~{fasta} ./fasta
        fasta=`ls ./fasta/*.fa`
        cp ~{fai} ./fasta/
        fai=`ls ./fasta/*.fai`
        cp ~{dict} ./fasta/
        dict=`ls ./fasta/*.dict`
        mkdir -p dbsnp
        cp ~{dbsnp} ./dbsnp/
        dbsnp=`ls ./dbsnp/*.vcf.gz`
        zcat "${dbsnp}" > ./${dbsnp%.gz}
        dbsnp=${dbsnp%.gz}
        gatk --java-options -Xmx40g \
            IndexFeatureFile \
            -I ${dbsnp%.gz}
        gatk --java-options -Xmx40g \
            GenotypeGVCFs \
                -R "${fasta}" \
                --dbsnp "${dbsnp}" \
                -V "~{gvcf}" \
                -new-qual -G StandardAnnotation \
                -O ~{sampleID}.vcf.gz \
                -isr INTERSECTION
    >>>
    runtime {
        dx_instance_type: "mem3_ssd1_v2_x8"
        docker: "quay.io/biocontainers/gatk4:4.2.0.0--0"
    }
    output {
        File vcf = "~{sampleID}.vcf.gz"
        File tbi = "~{sampleID}.vcf.gz.tbi"
    }
}