import pandas as pd
import argparse
import time
import datetime

t = datetime.datetime.now()
ts_string = t.strftime('%m%d%Y_%H%M%S')

parser = argparse.ArgumentParser()

parser.add_argument('--infile', '-i', type=str, help='input filename with coverage stats from VarSeq', required=True)
parser.add_argument('--outfile', '-o', type=str, help='output filename')
parser.add_argument('--mean', '-m', nargs='?', const='5.0', default='5.0', type=float, help='minimum average coverage for all exons, default = 5.0')
parser.add_argument('--min20', type=float, nargs='?', const='100.0', default='100.0', help='minimum value for 100 % 20x coverage, default = 100.0')
parser.add_argument('--min100', type=float, nargs='?', const='100.0', default='100.0', help='minimum value for 100% 100x coverage, default = 100.0')
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

min_mean_depth = float(args.mean)
min_20x_value = float(args.min20)
min_100x_value = float(args.min100)

df = pd.read_csv(args.infile, sep='\t')


if not args.geneList:
    print('no gene list supplied, using names derived from file.')
    nameCol = df['Name'].tolist()
    splitList = [x.split('/') for x in nameCol]
    nameList = list(set([x[0] for x in splitList]))
    geneList = nameList
    # print('using geneList:')
    # [print(x) for x in geneList]
else:
    geneList = []
    with open(args.geneList) as fOpen:
        for i in fOpen:
            i = i.rstrip('\r\n')
            geneList.append(i)

def getUniqueExons(df):
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
        return (False, [str(x[0]) + ':' + str(x[2]) for x in scannedList])
    else:
        # scannedValues = [x[2] for x in l]
        # update 10282021 JRC: include all exons with sufficient coverage
        # replaced code below that took max and returned one value 
        # with returning entire list of exons
        return (True, [str(x[0]) + ':' + str(x[2]) for x in l])



# columns for output DF
gene_name_list = []
mean_depth_bool_list = []
mean_depth_list = []
cov_20x_bool_list = []
cov_20x_list = []
cov_100x_bool_list = []
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
    for r in exonList:
        exonSublist = tmpdf.loc[tmpdf['Name'].str.endswith(r),:]
        # check mean depth
        meanDepthList = exonSublist[meanDepth_col].tolist()
        if max(meanDepthList) < min_mean_depth:
            meanDepth_covered_tuple.append((r, False, max(meanDepthList)))
        else:
            meanDepth_covered_tuple.append((r, True, max(meanDepthList)))
        # check min %20x coverage
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
    exonsWithLowMean = getCoverageStringAsList(meanDepth_covered_tuple)
    if exonsWithLowMean[0] == True:
        mean_depth_bool_list.append(True)
        outList = [str(x) for x in exonsWithLowMean[1]]
        outString = ','.join(outList)
        mean_depth_list.append(outString)
    else:
        mean_depth_bool_list.append(False)
        outList = [str(x) for x in exonsWithLowMean[1]]
        outString = ','.join(outList)
        mean_depth_list.append(outString)
        outString = 'low coverage: ' + outString
    # 20x results
    exonsWithLow20x = getCoverageStringAsList(min20x_covered_tuple)
    if exonsWithLow20x[0] == True:
        cov_20x_bool_list.append(True)
        outList = [str(x) for x in exonsWithLow20x[1]]
        outString = ','.join(outList)
        cov_20x_list.append(outString)
    else:
        cov_20x_bool_list.append(False)
        outList = [str(x) for x in exonsWithLow20x[1]]
        outString = ','.join(outList)
        cov_20x_list.append(outString)
    # 100x results
    exonsWithLow100x = getCoverageStringAsList(min100x_covered_tuple)
    if exonsWithLow100x[0] == True:
        cov_100x_bool_list.append(True)
        outList = [str(x) for x in exonsWithLow100x[1]]
        outString = ','.join(outList)
        cov_100x_list.append(outString)
    else:
        cov_100x_bool_list.append(False)
        outList = [str(x) for x in exonsWithLow100x[1]]
        outString = ','.join(outList)
        cov_100x_list.append(outString)
df_out = pd.DataFrame({
    'Gene_Name' : gene_name_list,
    # update 10282021 JRC: Boolean no longer needed 
    # 'Exons_Above_MeanCov_Threshold' : mean_depth_bool_list, 
    'Mean_Coverage_for_Exons' : mean_depth_list, 
    'Exons_at_100_20x' : cov_20x_bool_list, 
    'Exons_20x_Levels' : cov_20x_list,
    'Exons_at_100_100x' : cov_100x_bool_list,
    'Exons_100x_Levels' : cov_100x_list
})

df_out.to_csv(outFileName, sep='\t', index=False)
print(f'job done, created output file named {outFileName}')
