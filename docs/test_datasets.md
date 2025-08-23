# Test Datasets for Validating ConTra Analysis Pipeline

## Overview

This document identifies publicly available multi-omics datasets from model organisms with well-characterized epigenetic regulatory interactions. These datasets can be used to validate the ConTra analysis pipeline by comparing computed results against known biological relationships.

## Dataset Requirements

For effective validation, test datasets should include:

1. **Multi-omics data types** matching ConTra's inputs:
   - Gene expression (RNA-seq)
   - Long non-coding RNA (lncRNA) expression
   - MicroRNA (miRNA) expression
   - DNA methylation (WGBS/bisulfite-seq)

2. **Time-series or condition-dependent data** to enable context-dependent analysis

3. **Well-documented regulatory interactions** from literature

4. **Model organism** with extensive regulatory annotations

## Recommended Test Datasets

### 1. Arabidopsis thaliana - ENCODE Plant Project

**Dataset**: AtGenExpress + EpiGenome Project
- **Organism**: *Arabidopsis thaliana* (model plant)
- **Data types**: RNA-seq, lncRNA-seq, small RNA-seq, bisulfite-seq
- **Conditions**: Multiple developmental stages and stress conditions
- **Sample size**: ~50-100 samples across conditions
- **Time points**: Development stages (seedling, rosette, flowering)

**Known Interactions**:
- **FLC-miR172-FT flowering pathway**: FLC (flowering locus C) is regulated by DNA methylation and interacts with miR172 to control FT (flowering time)
- **Cold response**: CBF transcription factors are regulated by DNA methylation and coordinate with cold-responsive miRNAs
- **Circadian clock**: TOC1/CCA1 network involves lncRNAs and methylation-dependent regulation

**Access**: NCBI GEO, TAIR database
**Validation Value**: High - extensive literature on flowering time and stress response networks

### 2. Caenorhabditis elegans - modENCODE Project

**Dataset**: modENCODE developmental time series
- **Organism**: *C. elegans* (model nematode)  
- **Data types**: RNA-seq, small RNA-seq, ChIP-seq, DNA methylation
- **Conditions**: Developmental time course (embryo to adult)
- **Sample size**: ~40 samples across 4 developmental stages
- **Time points**: L1, L2, L3, L4 larval stages

**Known Interactions**:
- **let-7 miRNA pathway**: let-7 regulates developmental timing through lin-41, hbl-1 targets
- **lin-4 miRNA cascade**: Controls larval development timing through lin-14, lin-28
- **Germline silencing**: piRNAs and DNA methylation coordinate to silence repetitive elements

**Access**: ENCODE portal, WormBase
**Validation Value**: Very High - extremely well-characterized developmental network

### 3. Drosophila melanogaster - modENCODE Project

**Dataset**: modENCODE developmental time series
- **Organism**: *D. melanogaster* (model fly)
- **Data types**: RNA-seq, small RNA-seq, ChIP-seq, DNA methylation
- **Conditions**: Developmental stages (embryo, larva, pupa, adult)
- **Sample size**: ~30 samples across development
- **Time points**: 12 developmental time points

**Known Interactions**:
- **Hox gene regulation**: DNA methylation and miRNAs coordinate Hox gene expression
- **Wing development**: dpp/BMP pathway involves lncRNAs and epigenetic regulation
- **Sex determination**: Sex-lethal pathway includes miRNA and methylation components

**Access**: ENCODE portal, FlyBase
**Validation Value**: High - well-studied developmental networks

### 4. Mouse (Mus musculus) - ENCODE Project

**Dataset**: ENCODE mouse tissue atlas
- **Organism**: *M. musculus* (model mammal)
- **Data types**: RNA-seq, lncRNA-seq, miRNA-seq, WGBS
- **Conditions**: Multiple tissues and developmental stages  
- **Sample size**: ~80 samples across tissues/stages
- **Time points**: E11.5, E14.5, E16.5, P0 embryonic stages

**Known Interactions**:
- **X-chromosome inactivation**: Xist lncRNA, miRNAs, and DNA methylation
- **Imprinting clusters**: H19/Igf2 locus with DNA methylation and miRNA regulation
- **Pluripotency network**: Oct4, Sox2, Nanog with regulatory miRNAs and methylation

**Access**: ENCODE portal, MGI database
**Validation Value**: Very High - extensive regulatory network knowledge

## Recommended Starting Dataset: C. elegans modENCODE

**Rationale for C. elegans**:
1. **Manageable size**: Smaller genome reduces computational requirements
2. **Exceptional annotation**: Best-characterized regulatory networks among all organisms
3. **Clear temporal patterns**: Well-defined developmental stages with known transitions
4. **Established interactions**: Extensively validated miRNA and epigenetic regulatory cascades
5. **Data availability**: High-quality, standardized datasets available

## Expected Validation Results

For C. elegans developmental time series, ConTra should identify:

1. **let-7 pathway interactions**: 
   - let-7 miRNA negatively correlates with lin-41, hbl-1 expression
   - Context-dependent: strongest in L3-L4 transition
   - Methylation involvement in let-7 promoter regulation

2. **lin-4 temporal cascade**:
   - lin-4 miRNA anticorrelates with lin-14, lin-28
   - Context-dependent: L1-L2 transition specificity
   - lncRNA involvement in lin-4 regulation

3. **Germline regulatory network**:
   - piRNA abundance correlates with repetitive element silencing
   - DNA methylation marks associate with germline gene expression
   - Context-dependent: germline vs somatic tissue specificity

## Implementation Plan

1. **Download C. elegans modENCODE data**
2. **Reformat to ConTra input format** (match current Acropora structure)  
3. **Run ConTra analysis pipeline**
4. **Compare results against known interactions**
5. **Generate validation report**

## Data Access Information

**C. elegans modENCODE datasets**:
- **Portal**: https://www.encodeproject.org/
- **Organism filter**: *Caenorhabditis elegans*
- **Assay types**: RNA-seq, small RNA-seq, bisulfite-seq
- **Accessions**: ENCSR000XXX series (specific IDs to be determined)
- **File formats**: FASTQ, BAM, bigWig, BED

**Processing requirements**:
- Raw data preprocessing (if needed)
- Expression quantification
- Format conversion to match ConTra CSV structure
- Sample metadata alignment

## Success Criteria

A successful validation should show:

1. **Sensitivity**: ConTra detects known regulatory interactions (>80% of literature-validated interactions)
2. **Specificity**: Limited false positives (<20% of top interactions are novel/unvalidated)
3. **Context-dependency**: Temporal specificity matches known developmental windows
4. **Statistical significance**: P-values and effect sizes are appropriate for validated interactions
