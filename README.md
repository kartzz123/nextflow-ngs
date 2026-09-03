Containerized NGS Variant Calling PipelineA fully containerized, end-to-end Next-Generation Sequencing (NGS) variant calling workflow built with Nextflow (DSL2), Docker, and standard bioinformatics tools.This repository packages the environment, sample dataset (GRCh38 Chromosome 17), and workflow logic into a single container image that can run reproducibly anywhere Docker is installed.Pipeline OverviewThe Nextflow workflow inside this container runs the following sequential steps:Quality Control & Trimming (fastp): Cleans raw paired-end FASTQ reads and generates HTML/JSON QC reports.Reference Indexing (BWA & Samtools): Indexes the reference FASTA genome (Genome.fa).Alignment (BWA-MEM & Samtools): Aligns trimmed reads to the indexed reference and outputs a sorted, indexed BAM file.Variant Calling (BCFTools): Calls raw genomic variants from the alignment using mpileup and call.Normalization (BCFTools): Normalizes indels and left-aligns variants against the reference genome.Sorting & Indexing (BCFTools): Sorts the normalized VCF file and generates Tabix/CSI indices.Filtering (BCFTools): Filters low-quality variants (retaining variants with QUAL >= 20 and DP >= 10).Consensus Generation (BCFTools): Applies filtered variants to the reference genome to produce a consensus sequence FASTA file.PrerequisitesDocker: Installed and running on your system.Storage: At least 3 GB of free space (for the image build and GRCh38 Chromosome 17 reference download).Quick Start1. Build the Docker ImageBuild the image locally. This step automatically downloads the GRCh38 Chromosome 17 reference sequence and creates sample FASTQ input files.Bashdocker build -t ngs-pipeline:latest .
2. Run the ContainerRun the container to process the default sample files packaged inside the image:Bashdocker run --rm ngs-pipeline:latest
Running with Custom Host DataYou can pass your own reference genome and paired-end FASTQ files at runtime by mounting a local directory into the /data folder of the container.Bashdocker run --rm \
  -v /path/to/your/local_data:/data \
  -e READ1="/data/sample_R1.fastq.gz" \
  -e READ2="/data/sample_R2.fastq.gz" \
  -e REFERENCE="/data/my_reference.fa" \
  -e OUTDIR="/data/results" \
  ngs-pipeline:latest
Directory & Output StructureWhen the pipeline finishes running, output files are organized under the designated output directory (/data/results by default):Plaintextresults/
├── trimmed/          # Cleaned FASTQ files and fastp QC reports (.html/.json)
├── reference/        # Indexed reference FASTA files
├── aligned/          # Sorted BAM alignment files and .bai indices
├── vcf/
│   ├── raw/          # Unfiltered raw VCF outputs
│   ├── normalized/   # Left-aligned and normalized VCF outputs
│   ├── sorted/       # Sorted VCF files
│   └── filtered/     # Quality-filtered VCF files and Tabix indices (.tbi)
└── consensus/        # Final consensus sequence FASTA file (.consensus.fa)
To extract results directly to your local host machine, mount a local output path to /data/results:Bashdocker run --rm -v $(pwd)/my_results:/data/results ngs-pipeline:latest
Environment VariablesThe docker entrypoint script (run-ngs) supports the following configurable runtime environment variables:VariableDefault ValueDescriptionREAD1/data/ERR2356709_1.fastqPath to Forward Read 1 FASTQREAD2/data/ERR2356709_2.fastqPath to Reverse Read 2 FASTQREFERENCE/data/Genome.faPath to Reference FASTA GenomeOUTDIR/data/resultsPath to save final pipeline outputsToolchain Version InfoBase Image: Ubuntu 22.04 LTSWorkflow Engine: Nextflow (DSL2)Java Runtime: OpenJDK 17Bioinformatics Suite:fastpbwasamtoolsbcftools
