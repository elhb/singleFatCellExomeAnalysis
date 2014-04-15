
#
# Imports
#
import sys
import os

#
# get sample names
#
path = sys.argv[1]
walker = os.walk(path)
samples = walker.next()[1]

#
# Get the data and output table for MS Excel
#
data = { sample:{} for sample in samples }
outfile = open(path+'/mappingStats.tsv','w')
outfile.write('sample\ttotalReads\tmappingRate\tonTargetPercentage\tmeanTargetCoverage\tpercentageDuplication\n') # header for tsv file
for sample in samples:

    #
    # Read infiles
    #
    markDups = open(path+sample+'/'+sample+'.MarkDupsMetrix','r')
    hsMetrix = open(path+sample+'/'+sample+'.recalibrated.final.bam.hs_metrics.summary.txt','r')
    
    dupsLine = markDups.read().split('\n')[7].rstrip()
    hsLine = hsMetrix.read().split('\n')[7].rstrip()
    
    percentageDuplication   = float(dupsLine.split('\t')[7].replace(',','.'))
    totalReads              = int(hsLine.split('\t')[5])
    mappingRate             = float(hsLine.split('\t')[11].replace(',','.'))
    onTargetPercentage      = float(hsLine.split('\t')[17].replace(',','.'))
    meanTargetCoverage      = float(hsLine.split('\t')[21].replace(',','.'))
    
    markDups.close()
    hsMetrix.close
    
    #
    # save data to dictionary
    #
    data[sample]['percentageDuplication'] = percentageDuplication
    data[sample]['totalReads'] = totalReads
    data[sample]['mappingRate'] = mappingRate
    data[sample]['onTargetPercentage'] = onTargetPercentage
    data[sample]['meanTargetCoverage'] = meanTargetCoverage
    
    #
    # write to outfile
    #
    outfile.write( str(sample)+'\t'+str(totalReads).replace('.',',')+'\t'+str(mappingRate).replace('.',',')+'\t'+str(onTargetPercentage).replace('.',',')+'\t'+str(meanTargetCoverage).replace('.',',')+'\t'+str(percentageDuplication).replace('.',',')+'\n')

outfile.close() # close outfile

#
# Make graphics
#
