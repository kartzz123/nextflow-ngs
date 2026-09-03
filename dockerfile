FROM ubuntu:22.04

# Set non-interactive environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt \
    REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    NXF_HOME=/opt/nextflow

# Install System + Bioinformatics Tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    openssl \
    curl \
    wget \
    unzip \
    git \
    procps \
    gzip \
    bzip2 \
    xz-utils \
    python3 \
    openjdk-17-jre \
    bwa \
    samtools \
    bcftools \
    fastp \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create Required Pipeline Directories
RUN mkdir -p /pipeline /pipeline/work /data /data/results /opt/nextflow
WORKDIR /data

# Download GRCh38 Chromosome 17 and save as /data/Genome.fa
RUN curl -fsSL https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz | gunzip -c > /data/Genome.fa

# Create sample input FASTQ files for ERR2356709
RUN echo "@ERR2356709.1\nAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" > /data/ERR2356709_1.fastq && \
    echo "@ERR2356709.1\nTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" > /data/ERR2356709_2.fastq

# Write Nextflow Workflow Script
WORKDIR /pipeline
RUN cat > /pipeline/main.nf <<'NEXTFLOW'
nextflow.enable.dsl=2

params.read1     = "/data/ERR2356709_1.fastq"
params.read2     = "/data/ERR2356709_2.fastq"
params.reference = "/data/Genome.fa"
params.outdir    = "/data/results"

process TRIM {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: reads
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    """
    set -euo pipefail
    fastp -i ${r1} -I ${r2} \
        -o ${sample_id}_1.trimmed.fastq.gz -O ${sample_id}_2.trimmed.fastq.gz \
        --detect_adapter_for_pe --thread 4 \
        --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json
    """
}

process INDEX_REFERENCE {
    tag "Indexing Reference"
    publishDir "${params.outdir}/reference", mode: 'copy', overwrite: true

    input:
    path ref

    output:
    tuple path(ref), path("ref_indices/*"), emit: indexed_reference

    script:
    """
    set -euo pipefail
    mkdir -p ref_indices
    bwa index -p ref_indices/${ref.name} ${ref}
    samtools faidx ${ref}
    cp ${ref}.fai ref_indices/
    """
}

process ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy', overwrite: true

    input:
    tuple path(ref), path("ref_indices/*")
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    set -euo pipefail
    bwa mem -t 4 ref_indices/${ref.name} ${r1} ${r2} | samtools sort -@ 4 -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """
}

process VARIANT_CALLING {
    tag "${sample_id}"
    publishDir "${params.outdir}/vcf/raw", mode: 'copy', overwrite: true

    input:
    tuple path(ref), path("ref_indices/*")
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi"), emit: raw_vcf

    script:
    """
    set -euo pipefail
    bcftools mpileup -Ou -f ${ref} ${bam} | bcftools call -mv -Oz -o ${sample_id}.vcf.gz
    bcftools index -f ${sample_id}.vcf.gz
    """
}

process NORMALIZE_VCF {
    tag "${sample_id}"
    publishDir "${params.outdir}/vcf/normalized", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)
    path ref, stageAs: 'ref.fa'

    output:
    tuple val(sample_id), path("${sample_id}.normalized.vcf.gz"), path("${sample_id}.normalized.vcf.gz.csi"), emit: normalized_vcf

    script:
    """
    set -euo pipefail
    bcftools norm -f ${ref} -Oz -o ${sample_id}.normalized.vcf.gz ${vcf}
    bcftools index -f ${sample_id}.normalized.vcf.gz
    """
}

process SORT_VCF {
    tag "${sample_id}"
    publishDir "${params.outdir}/vcf/sorted", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.vcf.gz"), path("${sample_id}.sorted.vcf.gz.csi"), emit: sorted_vcf

    script:
    """
    set -euo pipefail
    bcftools sort ${vcf} -Oz -o ${sample_id}.sorted.vcf.gz
    bcftools index -f ${sample_id}.sorted.vcf.gz
    """
}

process FILTER_VCF {
    tag "${sample_id}"
    publishDir "${params.outdir}/vcf/filtered", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), path("${sample_id}.filtered.vcf.gz.tbi"), emit: filtered_vcf

    script:
    """
    set -euo pipefail
    bcftools view -i 'QUAL>=20 && INFO/DP>=10' -Oz -o ${sample_id}.filtered.vcf.gz ${vcf}
    bcftools index -t -f ${sample_id}.filtered.vcf.gz
    """
}

process CONSENSUS {
    tag "${sample_id}"
    publishDir "${params.outdir}/consensus", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)
    path ref, stageAs: 'ref.fa'

    output:
    path "${sample_id}.consensus.fa"

    script:
    """
    set -euo pipefail
    bcftools consensus -f ${ref} ${vcf} > ${sample_id}.consensus.fa
    """
}

workflow {
    reads_ch = channel.of(tuple("ERR2356709", file(params.read1, checkIfExists: true), file(params.read2, checkIfExists: true)))
    reference_ch = channel.of(file(params.reference, checkIfExists: true))

    trimmed_ch = TRIM(reads_ch)
    indexed_reference_ch = INDEX_REFERENCE(reference_ch)
    aligned_ch = ALIGN(indexed_reference_ch.indexed_reference, trimmed_ch.reads)
    raw_vcf_ch = VARIANT_CALLING(indexed_reference_ch.indexed_reference, aligned_ch.bam)
    normalized_vcf_ch = NORMALIZE_VCF(raw_vcf_ch.raw_vcf, reference_ch)
    sorted_vcf_ch = SORT_VCF(normalized_vcf_ch.normalized_vcf)
    filtered_vcf_ch = FILTER_VCF(sorted_vcf_ch.sorted_vcf)
    CONSENSUS(filtered_vcf_ch.filtered_vcf, reference_ch)
}
NEXTFLOW

# Create Entrypoint Execution Script
RUN cat > /usr/local/bin/run-ngs <<'SCRIPT'
#!/bin/bash
set -euo pipefail

R1="${READ1:-/data/ERR2356709_1.fastq}"
R2="${READ2:-/data/ERR2356709_2.fastq}"
REF="${REFERENCE:-/data/Genome.fa}"
OUT="${OUTDIR:-/data/results}"

if [ ! -f "$R1" ] || [ ! -f "$R2" ] || [ ! -f "$REF" ]; then
    echo "ERROR: Missing required input files."
    echo "Expected: $R1, $R2, $REF"
    exit 1
fi

mkdir -p "$OUT"

/usr/local/bin/nextflow run /pipeline/main.nf \
    --read1 "$R1" \
    --read2 "$R2" \
    --reference "$REF" \
    --outdir "$OUT" \
    -work-dir /pipeline/work
SCRIPT

RUN chmod +x /usr/local/bin/run-ngs
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/run-ngs"]
