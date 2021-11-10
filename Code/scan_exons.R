
library(tidyverse)
library(stringr)
library(gtools)

# required variables
##
inFile = c('r_file.sb.txt') # input file name
outFileName = c('r_file_sb.output.txt') # output file name

# functions
##
change_colnames <- function() {
  col1 <- c('Region')
  col2 <- c('Name')
  col3 <- paste0('Counted.Bases')
  col4 <- paste0('Mean.Depth')
  col5 <- paste0('MinDepth')
  col6 <- paste0('MaxDepth')
  col7 <- paste0('CovPercent.1x')
  col8 <- paste0('CovPercent.20x')
  col9 <- paste0('CovPercent.100x')
  col10 <- paste0('CovPercent.500x')
  col11 <- paste0('Mean.MQ')
  col12 <- paste0('Reads.in.Region')
  column_names = c(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12)
  return(column_names)
}

filterCoverage <- function(d_df, x_col) {
  filtered_df = d_df[which(d_df[, x_col] != 100.0),]
  filtered_df = filtered_df %>% mutate(cov_values = 100 - filtered_df[,x_col])
  exon_names = as.vector(filtered_df$ename)
  cov_values_vec = as.vector(filtered_df$cov_values)
  if (length(cov_values_vec) == 0) {
    return(c(''))
  } else {
    cov_values_vec = round(cov_values_vec, 3)
    cov_vec = c(sapply(seq_along(exon_names), function(i) paste0(exon_names[i], ":", cov_values_vec[i])))
    cov_vec_string = paste0(as.character(cov_vec), sep=",", collapse="")
    return(cov_vec_string)  
  }
}

# Workflow
##
df = read.csv(inFile, sep = '\t')

column_names = change_colnames()
colnames(df) <- column_names
df_data = df %>% separate(Name, sep='/', c("gname", "tname", "ename"))


# get unique list of gene names
gname_unique <- as.character(unique(df_data$gname))

# create output dataframe
df_output = setNames(
  data.frame(
    matrix(ncol=7, nrow=0)
    ),
  c("Gene_Name",
    "Transcript_Number",
    "Number_of_Exons",
    "Mean_Coverage_all_Exons",
    "Mean_Coverage_for_Each_Exon",
    "Exons_with_nucleotides_lessthan_20x_asPercent",
    "Exons_with_nucleotides_lessthan_100x_asPercent"
    )
)
# parse data into output dataframe
for (g in gname_unique) {
  # create subset dataframe then sort by exon
  tmpDf = df_data[which(df_data$gname == g),]
  tmpDf = tmpDf[mixedorder(tmpDf$ename),]
  # output column 1
  gname_vec = as.character(tmpDf[1, 'gname'])
  # output column 2
  tname_vec = as.character(tmpDf[1, 'tname'])
  # output column 3
  m_names = as.vector(tmpDf$ename)
  exonCt = length(m_names)
  # output column 4: mean coverage for all exons
  m_vec_float = mean(tmpDf$Mean.Depth)
  # output column 5: mean coverage column each exon
  m_vals = as.vector(tmpDf$Mean.Depth)
  m_vec = c(sapply(seq_along(m_vals), function(i) paste0(m_names[i], ":", m_vals[i])))
  m_vec_string = paste0(as.character(m_vec), sep=",", collapse="")
  # output column 6: low covergae at 20x
  covCheck_20x = filterCoverage(tmpDf, 'CovPercent.20x')
  # output column 7: low coverage at 100x
  covCheck_100x = filterCoverage(tmpDf, 'CovPercent.100x')
  # add to output dataframe
  t_df_output = data.frame(
    Gene_Name = gname_vec,
    Transcript_Number = tname_vec,
    Number_of_Exons = exonCt,
    Mean_Coverage_all_Exons = m_vec_float,
    Mean_Coverage_for_Each_Exon = m_vec_string,
    Exons_with_nucleotides_lessthan_20x_asPercent = covCheck_20x,
    Exons_with_nucleotides_lessthan_100x_asPercent = covCheck_100x
  )
  df_output = bind_rows(df_output, t_df_output)
}
# Exons_with_nucleotides_lessthan_20x (:% of that exon covered to 20x or less (1- numbers)
final_colNames = c("Gene_Name",
                   "Transcript_Number",
                   "Number_of_Exons",
                   "Mean_Coverage_all_Exons",
                   "Mean_Coverage_for_Each_Exon",
                   "Exons_with_nucleotides_lessthan_20x_asPercent_of_exon_covered_to_20x_or_less",
                   "Exons_with_nucleotides_lessthan_100x_asPercent_of_exon_covered_to_100x_or_less"
)
colnames(df_output) <- final_colNames

# output to tab-delimited file
write.table(df_output, file=outFileName, sep = "\t")
outputDone = paste0('created output file ', outFileName)
print(outputDone)