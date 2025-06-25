# PMS2-Variant-Analysis-Pipeline
Dockerized implementation of PMS2_vaR pipeline for PMS2 variant analysis.  Based on Munté et al. (2024) - enhanced with containerization and improved result compilation.

## Overview

This pipeline addresses the challenge of distinguishing between PMS2 gene variants and its pseudogene (PMS2CL) variants, which share high sequence homology (>98% identity). The pipeline uses a sophisticated two-approach methodology to improve variant calling accuracy while minimizing the need for labor-intensive long-range PCR validation.

## Original Publication

Based on the work by Munté et al. (2024):
- **Paper**: "Open-Source Bioinformatic Pipeline to Improve PMS2 Genetic Testing Using Short-Read NGS Data"
- **Journal**: Journal of Molecular Diagnostics, Vol. 26, No. 8, August 2024
- **DOI**: https://doi.org/10.1016/j.jmoldx.2024.05.005
- **Original Repository**: https://github.com/emunte/PMS2_vaR

## Key Features

- **Two-Approach Methodology**: Combines general approach (masked reference alignment) with E11-specific approach (invariant position filtering)
- **Automated Decision Algorithm**: Provides recommendations for which variants need LR-PCR validation
- **Dockerized Environment**: All dependencies pre-installed and configured
- **Enhanced Result Compilation**: Comprehensive analysis and reporting of variant calls

## Requirements

- Docker
- FASTQ files (paired-end sequencing data)
- ~24GB RAM recommended

## Quick Start

### 1. Build Docker Image

```bash
docker build -t pms2_pipeline .
```

### 2. Prepare Data Structure

Create the following directory structure with your FASTQ files:
```
project_directory/
├── config/
│   ├── tools.yaml
│   └── vardicjavaParams.yaml
├── data/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── sample2_R2.fastq.gz
├── PMS2_vaR/
├── reference/
└── results/
```

### 3. Run Pipeline

```bash
# Start Docker container
docker run -v $(pwd):/home/rstudio/data --rm -it -m 24g --memory-swap 24g pms2_pipeline bash

# Create original BAM files from FASTQ
for R1_FILE in data/*_R1.fastq.gz; do
    R2_FILE="${R1_FILE/_R1/_R2}"
    SAMPLE_NAME=$(basename ${R1_FILE/_R1.fastq.gz/})
    
    bwa mem -t 4 -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA" \
        /home/rstudio/reference/chr7.fa \
        ${R1_FILE} ${R2_FILE} | \
        samtools view -bS - | \
        samtools sort -@ 4 -o data/${SAMPLE_NAME}.original.bam
    
    samtools index data/${SAMPLE_NAME}.original.bam
done

# Create BAM file list
ls -1 $PWD/data/*.original.bam > data/bam_files.txt

# Run PMS2_vaR pipeline
cd /home/rstudio/PMS2_vaR
Rscript run_PMS2_vaR.R \
  -t /home/rstudio/data/config/tools.yaml \
  -b /home/rstudio/data/data/bam_files.txt \
  -r /home/rstudio/reference/modified_reference/chr7_masked.fa \
  -g hg38 \
  -v /home/rstudio/data/config/vardicjavaParams.yaml \
  -o /home/rstudio/data/results \
  -n PMS2_analysis

# Compile and analyze results
cd /home/rstudio/data
Rscript compile_results_enhanced.R results/PMS2_analysis_$(date +%Y-%m-%d)
```

## Configuration Files

### tools.yaml
```yaml
samtools: /usr/bin/samtools
picard: /opt/tools/picard.jar
bwa: /usr/bin/bwa
vardict: /usr/local/bin/VarDict
```

### vardicjavaParams.yaml
```yaml
c: 1
S: 2
E: 3
g: 4
f: 0.01
r: 4
q: 20
Q: 10
V: 50
m: 8
th: 10
```

## Output Files

- **Individual Results**: `results/PMS2_analysis_DATE/final_results/SAMPLE_PMS2_variants.txt`
- **Summary Report**: `results/PMS2_analysis_DATE/enhanced_variant_summary.csv`
- **Pathogenic Variants**: `results/PMS2_analysis_DATE/pathogenic_variants.csv` (if any found)

## Pipeline Methodology

1. **Reference Preparation**: Creates masked reference genome with PMS2CL region replaced by N's
2. **BAM Preprocessing**: Aligns FASTQ files to unmasked reference genome
3. **Two-Approach Analysis**:
   - **General Approach**: Forces alignment of PMS2/PMS2CL reads to PMS2 reference
   - **E11 Approach**: For exon 11, uses only reads overlapping invariant paralogous sequence variants
4. **Variant Calling**: Uses VarDictJava for variant detection
5. **Classification**: Applies decision algorithm to recommend LR-PCR validation

## Key Differences from Original

- **Dockerized Environment**: Complete containerization for reproducibility
- **Enhanced Result Compilation**: Additional analysis script for comprehensive reporting
- **Streamlined Workflow**: Optimized for FASTQ input files
- **Updated Dependencies**: Compatible with latest tool versions

## Citation

If you use this pipeline in your research, please cite:

Munté, E., Feliubadaló, L., Del Valle, J., et al. Open-Source Bioinformatic Pipeline to Improve PMS2 Genetic Testing Using Short-Read NGS Data. *J Mol Diagn* 2024, 26: 727-738. https://doi.org/10.1016/j.jmoldx.2024.05.005

## Support

For issues related to this dockerized implementation, please open an issue in this repository.
For questions about the original pipeline methodology, refer to the original publication and repository.
