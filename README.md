
# ATACseq_Pipeline

This directory contains a Snakemake workflow for processing bulk ATAC-seq data. The pipeline automates quality control, trimming, alignment, filtering, peak calling, and summarization for multiple samples.

## Workflow Steps

1. **FastQC**: Quality control of raw FASTQ files.
2. **Trimmomatic**: Adapter and quality trimming of reads.
3. **Bowtie2**: Alignment of trimmed reads to a reference genome.
4. **Samtools**: Sorting, filtering, and indexing of BAM files.
5. **MACS2**: Peak calling for accessible chromatin regions.
6. **MultiQC**: Aggregated report of QC and quantification results.

## Directory Structure
- `Snakefile`: Main workflow definition.
- `config.yaml`: Configuration file with paths and parameters.
- `ARMOR/`: Contains the ARMOR pipeline for RNA-seq analysis (optional, see `ARMOR/README.md`).
- `bin/`: Executable scripts and binaries used in the workflow.
- `data/`: Raw FASTQ files and related data.
- `logs/`: Log files for each step.
- `output/`: Processed data outputs (alignments, peaks, QC reports).
- `proj_src/`: Project-specific source code and scripts.
- `results/`: Final results and analysis outputs.
- `src/`: Additional scripts and source files for data processing and analysis.
- `tmp/`: Temporary files generated during workflow execution.

## Usage

### 1. Prerequisites

Make sure you have Snakemake installed. You can install it using conda:

```bash
conda install -c conda-forge -c bioconda snakemake
```

### 2. Configuration
Edit `config.yaml` to set paths and parameters for your data and references.
- Sample names
- Input/output paths
- Reference file locations
- Tool parameters

### 3. Running the workflow

#### Option A: Submit to SLURM cluster
```bash
sbatch submit_snakemake.sh
```

#### Option B: Run locally (for testing)
```bash
snakemake --cores 8 --use-conda
```

#### Option C: Dry run (to check workflow)
```bash
snakemake --dry-run
```

### 4. Workflow visualization

Generate a workflow diagram:
```bash
snakemake --dag | dot -Tpng > workflow.png
```

### 5. Output
1. **FastQC**: Quality control reports in `output/fastqc/`.
2. **Results**: Outputs will be found in the `output/` and `results/` directories. Main peak calling results are in `output/macs2/`.

## Requirements
- Snakemake
- Modules: fastqc, trimmomatic, bowtie2, samtools, macs2, multiqc, singularity
- Cluster environment (recommended)

## Customization
- Adjust sample detection, references, and tool parameters in `config.yaml`.
- Modify cluster submission scripts for resource allocation.

### For different library types:
- **Paired-end or single-end**: Change parameters in `config.yaml` accordingly.
- **Different reference genomes**: Modify reference paths in `config.yaml`.

## Key Features

1. **Modular design**: Each step is a separate rule
2. **Dependency management**: Snakemake automatically handles job dependencies
3. **Parallel execution**: Multiple samples can be processed simultaneously
4. **Configuration-driven**: Easy to modify parameters without editing the main workflow
5. **Resource management**: Better integration with SLURM scheduler
6. **Reproducibility**: Workflow tracks input/output dependencies

## Troubleshooting

1. Check SLURM job status: `squeue -u $USER`
2. View workflow status: `snakemake --summary`
3. Check individual rule logs in the SLURM output files and `logs/` directory

# TODO
1. Summarize peak counts and generate a final report
2. Integrate downstream analysis (e.g., motif analysis, differential accessibility)

## Contact
For questions or issues, contact: kstachel@uci.edu or refer to documentation in the `doc/` folder.

---

*Last updated: August 29, 2025*
