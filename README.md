# EpiNano_Koculi_Processing

### Epinano Data Processing and Plotting Pipeline for E. Coli ribosomal RNA using EpiNano output data for Nanopore Direct RNA Sequencing.

## Overview
This repository contains a set of scripts to process Epinano data and generate plots for RNA modifications and gene differences in 16s and 23s rRNA.

## Scripts

### 1. `rRNA_run.py`
This script serves as the main pipeline to process raw Epinano data and generate plots. It integrates the functionality of `process_epinano_replicates.py` and `plot_epinano_replicates.py`.

#### Usage
```bash

python rRNA_run.py --file1 <path_to_raw_data1> --file2 <path_to_raw_data2> [options]

```

## Arguments

    --file1: Path to the first raw data CSV file.
    --file2: Path to the second raw data CSV file.
    --file3: Filename for the combined output (default: output.csv).
    --processed_output: Filename for the processed output (default: processed_data.csv).
    --label1: Label for the first processed data (default: Replicate 1).
    --label2: Label for the second processed data (default: Replicate 2).
    --cap: Apply cap to the plots (optional).

### 2. process_epinano_replicates.py

This script processes and combines two Epinano CSV files, labels positions based on known modifications and gene differences, and saves the combined data.

```bash
python process_epinano_replicates.py --mode <combine|label> [options]
```

## Arguments

    --replicate1: Path to the first processed data CSV file.
    --replicate2: Path to the second processed data CSV file.
    --label1: Label for the first processed data (default: Replicate 1).
    --label2: Label for the second processed data (default: Replicate 2).
    --cap: Apply cap to the plots (optional).


## Example Workflow

1) Process raw data files:
```bash
python rRNA_run.py --file1 path/to/raw_data1.csv --file2 path/to/raw_data2.csv
```

2) Generate plots from processed data:

```bash
python rRNA_run.py --file1 path/to/raw_data1.csv --file2 path/to/raw_data2.csv --cap
```
