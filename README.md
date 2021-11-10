# Parse VarSeq Coverage Output
This repository has python and R code to convert VarSeq Coverage tab-delimited text output to a tab-delimited text file of coverage by exons.

# Usage
Both scripts require as input the text file from VarSeq, and both output a tab-delimited text file. The python script has added requirements and options, and is callable from commandline (Bash). The R script can be run via RStudio or commandline (Bash), and requires you to add the input file name and the output file name before the script will run.

## Python Usage

The [python script](https://github.com/disulfidebond/Parse_VarSeq_Coverage_Output/blob/main/Code/scan_exons.py) outputs low coverage exons from a provided input gene list and VarSeq Coverage text file via python `scan_exons.py`. The options are:

    --infile/-i the tab-delimited coverage file from VarSeq (required)
    --outfile the output file name
    --sampleName the sample name from the coverage file from VarSeq, must match exactly (required)
    --geneList the newline-delimited gene list file
    --force force the coverage file to run even if there is a mismatch in the provided gene names and the gene names in the coverage report
    
The input file is self explanatory. It is strongly recommended to provide an output file name, but the script will output the output file name to STDOUT.

The sampleName is the sample name from the coverage file from VarSeq. Special characters are allowed, but it *must* match exactly. For best results, copy and paste from the coverage file.

If no text gene list is provided, then the script will assume the gene names in the 'Names' column of the coverage file are correct, and will use those. If a text gene list is provided, then the script will verify the gene list against the gene names in the 'Names' column of the coverager file. If there is a mismatch between the provided gene list text file and the genes listed in the coverage file, the script will output a warning and stop.

If the `--force` option is used, then the script will output a warning if there is a mismatch between the provided gene list text file and the genes listed in the coverage file, but will still run. If no text gene list is provided, then this option is ignored.

## R Usage
The [R script]() can be run from RStudio or Bash. Before running, modify the `inFile` variable to hold the input file name, and the `outFileName` variable to hold the output file name. The script will not check if a file is present before overwriting it.

To run via commandline, you could do `Rscript scan_exons.py`, and to run it in RStudio, simply load the file and modify the variables as indicated.
