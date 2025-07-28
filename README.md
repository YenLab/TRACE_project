# **RUNX1 Utilizes Novel 3bp-Spaced Motif Arrays to Regulate Cell-Context-Dependent Function**

**:writing_hand:Author**: Shan Liu

**:email:Email**: 3shanliu3@gmail.com

Copyright (c) 2025 YenLab@SKLEH. All rights reserved.

## Introduction

The transcription factor RUNX1 orchestrates hematopoietic differentiation but drives leukemogenesis, creating a functional paradox. RUNX1 binds open chromatin via a canonical 12-bp motif but most in vivo targets lack this sequence, challenging conventional TF-DNA interaction models. 

To resolve this, we developed TRACE (Tagging and Recovery of Chromatin-associated Elements), a mammalian system employing endogenous DNA libraries to map RUNX1 binding in physiological contexts. 

We uncover: 1) Novel motif arrays: RUNX1 targets inaccessible chromatin through novel 3bp-spaced motif arrays, forming alternative binding platforms. 2) Epigenetic plasticity: DNA methylation at motif array regions inversely correlates with differential TET2 expression across HSPC and AML, revealing cell-type-specific epigenetic crosstalk. 3) Functional complexity: Motif arrays regulate both lineage-specific and leukemic programs, connecting architecture to pathological switching. 

Integrating TRACE-defined motif arrays with methylation dynamics, we demonstrate how RUNX1 toggles between hematopoiesis and leukemogenesis via a syntax-based mechanism, merging sequence arrangement with epigenetic context. 

<img align="left" width=450 src="https://github.com/SANGO63/TRACE_project/blob/main/graphic_summary/TRACE_Graphical_Summary.png">  

### :file_folder:Scripts organization

- Scripts head with `0_` is for general preprocessing NGS data, include

  `TRACE`, `ChIP-seq`,  `CUT&Tag`,  `WGBS`,  `MNase-seq` ,`RNA-seq` 

- Script head with series number is pipeline for specific analysis follows the order in the TRACE paper.

  Specifically, the `bash script` is for processing data and corresponding `R scripts` is used for downstream analysis, statistics and plotting.

- Published HT-SELEX data was processed using the Inomotif package with modified Python code, as detailed in the directory: Inomotif_modified
  
- Scripts in the Script_pipeline folder was called by pipelines.
