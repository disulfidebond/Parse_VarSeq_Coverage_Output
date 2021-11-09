# Parse VarSeq Coverage Output
This python script outputs low coverage exons from a provided input gene list and VarSeq Coverage text file.

usage is:

    --infile/-i the tab-delimited coverage file from VarSeq (required)
    --outfile the output file name
    --sampleName the sample name from the coverage file from VarSeq, must match exactly (required)
    --geneList the newline-delimited gene list file
    --force force the coverage file to run even if there is a mismatch in the provided gene names and the gene names in the coverage report
    
The input file is self explanatory. It is strongly recommended to provide an output file name, but the script will output the output file name to STDOUT.

The sampleName is the sample name from the coverage file from VarSeq. Special characters are allowed, but it *must* match exactly. For best results, copy and paste from the coverage file.

If no text gene list is provided, then the script will assume the gene names in the 'Names' column of the coverage file are correct, and will use those. If a text gene list is provided, then the script will verify the gene list against the gene names in the 'Names' column of the coverager file. If there is a mismatch between the provided gene list text file and the genes listed in the coverage file, the script will output a warning and stop.

If the `--force` option is used, then the script will output a warning if there is a mismatch between the provided gene list text file and the genes listed in the coverage file, but will still run. If no text gene list is provided, then this option is ignored.


# Integration with R

In addition, the code to run it in R is included.
