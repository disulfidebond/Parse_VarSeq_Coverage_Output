import pandas as pd
import argparse
import time
import datetime
import statistics
import sys

t = datetime.datetime.now()
ts_string = t.strftime('%m%d%Y_%H%M%S')

parser = argparse.ArgumentParser()

parser.add_argument('--infile', '-i', type=str, help='input filename with coverage stats from VarSeq', required=True)
parser.add_argument('--outfile', '-o', type=str, help='output filename')
parser.add_argument('--sampleName', '-s', type=str, help='sample name in VarSeq output file', required=True)
parser.add_argument('--geneList', '-g', type=str, help='newline-delimited text file with genes of interest, default is to infer gene names from the VarSeq file.')

args = parser.parse_args()
outFileName = ''
if not args.outfile:
    outFileName = 'scanned_output.' + ts_string + '.txt'
else:
    outFileName = str(args.outfile) + '.' + ts_string + '.txt'

geneList = ''
sampleName = args.sampleName

print(f'using sample name {sampleName}')

min_20x_value = 100.0
min_100x_value = 100.0

df = pd.read_csv(args.infile, sep='\t')


if not args.geneList:
    print('no gene list supplied, using names derived from file.')
    nameCol = df['Name'].tolist()
    splitList = [x.split('/') for x in nameCol]
    nameList = list(set([x[0] for x in splitList]))
    geneList = nameList
else:
    geneList = []
    with open(args.geneList) as fOpen:
        for i in fOpen:
            i = i.rstrip('\r\n')
            geneList.append(i)
    # update: compare provided gene list with list of genes, stop if mismatch
    check_geneList = df['Name'].tolist()
    checked_geneList = [x for x in check_geneList if x not in geneList]
    if len(checked_geneList) > 0:
        print('Warning, found genes in VarSeq that were not present in provided gene list:')
        print(checked_geneList)
        print('Exiting now.')
        sys.exit()

def getUniqueExons(df):
    # sort the exons in order by extracting the exon int, 
    # sorting, then adding back the prefix string
    # match is to the 'exonN' from the original dataframe
    # Dev Comment: this is kept in case having multiple exons and not canonical exons only is required
    geneList = df['Name'].tolist()
    exonList = [x.split('/') for x in geneList]
    exonList = [x[-1] for x in exonList]
    exonList = list(set(exonList))
    exonIdxList = [x.replace('ex', '') for x in exonList]
    exonIdxList = [int(x) for x in exonIdxList]
    exonIdxList.sort()
    exonList = [('ex' + str(x)) for x in exonIdxList]
    return exonList
def getCoverageStringAsList(l):
    scannedList = [x for x in l if x[1] == False]
    if len(scannedList) > 0:
        # update: get 20x || 100x coverage for exons 
        # that have less than 100% 20x || 100x coverage
        lowExonCoverageValues = [x[2] for x in scannedList]
        lowExonCoveragePercent = [(100-x) for x in lowExonCoverageValues]
        returnedList = []
        for idx in range(len(scannedList)):
            s = str(scannedList[idx][0]) + ':' + str(lowExonCoveragePercent[idx])
            returnedList.append(s)
        return (False, returnedList)
    else:
        return (True, [str(x[0]) + ':' + str(x[2]) for x in l])



# columns for output DF
gene_name_list = []
transcriptList = []
mean_depth_list = []
exon_ct_list = []
exons_mean_cov = []
cov_20x_list = []
cov_100x_list = []
# end columns definitions

# scan data, then format output
for n in geneList:
    gene_name_list.append(n)
    n_identifier = n + '/'
    tmpdf = df.loc[df['Name'].str.startswith(n_identifier),:]
    exonList = getUniqueExons(tmpdf)
    meanDepth_covered_tuple = []
    min20x_covered_tuple = []
    min100x_covered_tuple = []
    meanDepth_col = sampleName + ' Mean Depth'
    min20x_col = sampleName + ' % 20x'
    min100x_col = sampleName + ' % 100x'
    exonCt = 0
    for r in exonList:
        exonCt += 1
        exonSublist = tmpdf.loc[tmpdf['Name'].str.endswith(r),:]
        # get transcript name
        tNameList = tmpdf['Name'].tolist()
        tNameList = [x.split('/') for x in tNameList]
        transcriptName = tNameList[0][1]
        # get mean depth
        meanDepthList = exonSublist[meanDepth_col].tolist()
        meanDepth_covered_tuple.append((r, transcriptName, max(meanDepthList)))
        min20xList = exonSublist[min20x_col].tolist()
        if max(min20xList) < min_20x_value:
            min20x_covered_tuple.append((r, False, max(min20xList)))
        else:
            min20x_covered_tuple.append((r, True, max(min20xList)))
        # check min %100x coverage
        min100xList = exonSublist[min100x_col].tolist()
        if max(min100xList) < min_100x_value:
            min100x_covered_tuple.append((r, False, max(min100xList)))
        else:
            min100x_covered_tuple.append((r, True, max(min100xList)))
    # mean coverage results
    transcriptNameList = [x[1] for x in meanDepth_covered_tuple]
    transcriptNameList = list(set(transcriptNameList))
    transcriptName = transcriptNameList[0]
    if len(transcriptNameList) > 1:
        print(f'multiple transcripts detected for gene {n}:')
        print(transcriptNameList)
        transcriptName = '_'.join(transcriptNameList)
    # calculate arithm mean for all exons
    arithm_mean_list = [float(x[2]) for x in meanDepth_covered_tuple]
    arithm_mean = statistics.mean(arithm_mean_list)
    exonsWithMean = getCoverageStringAsList(meanDepth_covered_tuple)
    outList = [str(x) for x in exonsWithMean[1]]
    outString = ','.join(outList)
    # append to lists
    exons_mean_cov.append(arithm_mean)
    mean_depth_list.append(outString)
    transcriptList.append(transcriptName)
    exon_ct_list.append(exonCt)
    # 20x results
    exonsWithLow20x = getCoverageStringAsList(min20x_covered_tuple)
    if exonsWithLow20x[0] == True:
        outString = ''
        cov_20x_list.append(outString)
    else:
        outList = [str(x) for x in exonsWithLow20x[1]]
        outString = ','.join(outList)
        cov_20x_list.append(outString)
    # 100x results
    exonsWithLow100x = getCoverageStringAsList(min100x_covered_tuple)
    if exonsWithLow100x[0] == True:
        outString = '' 
        cov_100x_list.append(outString)
    else:
        outList = [str(x) for x in exonsWithLow100x[1]]
        outString = ','.join(outList)
        cov_100x_list.append(outString)
df_out = pd.DataFrame({
    'Gene_Name' : gene_name_list,
    'Transcript_Number' : transcriptList,
    'Number_of_Exons' : exon_ct_list, 
    'Mean_Coverage_All_Exons' : exons_mean_cov,
    'Mean_Coverage_for_Each_Exon' : mean_depth_list,  
    'Exons_with_nucleotides_lessthan_20x (:% of that exon covered to 20x or less (1- numbers)' : cov_20x_list,
    'Exons_with_nucleotides_lessthan_100x (:% of that exon covered to 100x or less (1- numbers)' : cov_100x_list
})

df_out.to_csv(outFileName, sep='\t', index=False)
print(f'job done, created output file named {outFileName}')
