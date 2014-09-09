

def main():

    #
    # Imports
    #
    import sys
    import os
    import re
    import gzip

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
    header = ['Sample','OriginalReads', 'Originalbases', 'TrimmedBasesWGAadapter','%WGAadapter', 'BasesCutadapt', 'TrimmedBasesCutadapt','%Cutadapt','Qualityprocessedbases','qualitytrimmedbases','%qualtrimmedbases','pairs','overallalignmentrateByBowtie','mappedpassingfilter','%mappedpassingfilter','properlypaired','%properlypaired','percentageDuplication','onTargetPercentage','meanTargetCoverage']
    print '\t'.join(header)
    outfile.write('\t'.join(header)+'\n')
#    outfile.write('sample\ttotalReadsPreFilter\ttotalReads\tmappingRate\tonTargetPercentage\tmeanTargetCoverage\tpercentageDuplication\n') # header for tsv file
    for sample in samples:
    
        print 'processing',sample,
        
        try:
            r1wgaTrimmingStats = open(path+sample+'/r1wgaTrimming.log.txt','r')
            line = r1wgaTrimmingStats.readline().rstrip()
            line = r1wgaTrimmingStats.readline().rstrip()
            r1Reads = int(line.split('\t')[1])
            line = r1wgaTrimmingStats.readline().rstrip()
            r1Bases = int(line.split('\t')[1])
            line = r1wgaTrimmingStats.readline().rstrip()
            r1TrimmedBases = int(line.split('\t')[1])
    
            r2wgaTrimmingStats = open(path+sample+'/r2wgaTrimming.log.txt','r')
            line = r2wgaTrimmingStats.readline().rstrip()
            line = r2wgaTrimmingStats.readline().rstrip()
            r2Reads = int(line.split('\t')[1])
            line = r2wgaTrimmingStats.readline().rstrip()
            r2Bases = int(line.split('\t')[1])
            line = r2wgaTrimmingStats.readline().rstrip()
            r2TrimmedBases = int(line.split('\t')[1])
            
            assert r1Reads == r2Reads,str(r1Reads)+' == '+str(r2Reads)
        #assert r1Reads == preTrimmingReadCount,str(r1Reads)+' == '+str(preTrimmingReadCount)
        #assert r1Bases+r2Bases==preTrimmingBaseCount, str(r1Bases)+'+'+str(r2Bases)+'=='+str(preTrimmingBaseCount)

        #print r1TrimmedBases,r2TrimmedBases, r1TrimmedBases+r2TrimmedBases,'bases trimmed by wga adapter trimming'
        except (IOError, IndexError):
            r1Bases = 'NA'
            r1Reads = 'NA'
            r2Bases = 'NA'
            r2Reads = 'NA'
            r1TrimmedBases = 'NA'
            r2TrimmedBases = 'NA'
            
            # get stats for original files
            
            trimmingScript = open(path+sample+'/script/'+sample+'.trimming.sh','r')
            while True:
                line = trimmingScript.readline().rstrip()
                if line == "": notFound=True;break
                if re.search('wgaAdapterTrimmer.py',line):
                    notFound=False
                    r1 = line.split(' ')[2]
                    line = trimmingScript.readline().rstrip()
                    r2 = line.split(' ')[2]
                    break
            trimmingScript.close()
            if notFound:
                trimmingScript = open(path+sample+'/script/'+sample+'.trimming.sh','r')
                while True:
                    line = trimmingScript.readline().rstrip()
                    if line == "": notFound=True;break
                    if re.search('cutadapt',line):
                        notFound=False
                        r1 = line.split()[11]
                        line = trimmingScript.readline().rstrip()
                        r2 = line.split()[11]
                        break                
                trimmingScript.close()
            if not notFound:
                r1 = gzip.open(r1,'r')
                r2 = gzip.open(r2,'r')
                linecount = [None,0,0]#[None,bufcount(r1.name),bufcount(r2.name)]
                r1.readline()
                r2.readline()
                seq1=r1.readline().rstrip()
                seq2=r2.readline().rstrip()
                assert len(seq1)==len(seq2)
                assert linecount[1]==linecount[2]
                preTrimmingReadCount = int(linecount[1])/4
                preTrimmingBaseCount = int(linecount[1])/4*int(len(seq1))*2
                #print preTrimmingBaseCount,'bases', preTrimmingReadCount,'reads',len(seq2),'read length'
                r1Reads = preTrimmingReadCount
                r2Reads = preTrimmingReadCount
                r1Bases = preTrimmingBaseCount/2
                r2Bases = preTrimmingBaseCount/2


        try:
            r1adapterTrimmingStats = open(path+sample+'/r1adapterTrimming.log.txt','r')
            for i in range(5):line = r1adapterTrimmingStats.readline().rstrip()
            r1ReadsCutadapt = int(line.split()[2].split(' ')[0])
            line = r1adapterTrimmingStats.readline().rstrip()
            r1BasesCutadapt = int(line.split()[2].split(' ')[0])
            line = r1adapterTrimmingStats.readline().rstrip()
            r1TrimmedReadsCutadapt = int(line.split()[2].split(' ')[0])
            line = r1adapterTrimmingStats.readline().rstrip()
            r1TrimmedBasesCutadapt = int(line.split()[2].split(' ')[0])
            
            r2adapterTrimmingStats = open(path+sample+'/r2adapterTrimming.log.txt','r')
            for i in range(5):line = r2adapterTrimmingStats.readline().rstrip()
            r2ReadsCutadapt = int(line.split()[2].split(' ')[0])
            line = r2adapterTrimmingStats.readline().rstrip()
            r2BasesCutadapt = int(line.split()[2].split(' ')[0])
            line = r2adapterTrimmingStats.readline().rstrip()
            r2TrimmedReadsCutadapt = int(line.split()[2].split(' ')[0])
            line = r2adapterTrimmingStats.readline().rstrip()
            r2TrimmedBasesCutadapt = int(line.split()[2].split(' ')[0])
    
            #assert preTrimmingBaseCount-r1TrimmedBases-r2TrimmedBases == r1BasesCutadapt+r2BasesCutadapt
            #print r1TrimmedBasesCutadapt + r2TrimmedBasesCutadapt,'bases trimmed out of', r1BasesCutadapt+r2BasesCutadapt,'bases'
        except  (IOError, IndexError):
            r1BasesCutadapt = 'NA'
            r1ReadsCutadapt = 'NA'
            r1TrimmedBasesCutadapt = 'NA'
            r1TrimmedReadsCutadapt = 'NA'
            r2BasesCutadapt = 'NA'
            r2ReadsCutadapt = 'NA'
            r2TrimmedBasesCutadapt = 'NA'
            r2TrimmedReadsCutadapt = 'NA'

        try:
            r1qualityTrimmingStats = open(path+sample+'/r1qualityTrimming.log.txt','r')
            line = r1qualityTrimmingStats.readline().rstrip()
            r1qualbases = int(line.split()[0])
            line = r1qualityTrimmingStats.readline().rstrip()
            r1qualtrimbases = int(line.split()[0])    
            r2qualityTrimmingStats = open(path+sample+'/r2qualityTrimming.log.txt','r')
            line = r2qualityTrimmingStats.readline().rstrip()
            r2qualbases = int(line.split()[0])
            line = r2qualityTrimmingStats.readline().rstrip()
            r2qualtrimbases = int(line.split()[0])
        except  (IOError, IndexError):
            r1qualbases = 'NA'
            r1qualtrimbases = 'NA'
            r1qualtrimbases = 'NA'
            r2qualbases = 'NA'

        try:
            removeEmptyreads =  open(path+sample+'/removeEmptyReads.log.txt','r')
            line = removeEmptyreads.readline().rstrip()
            line = removeEmptyreads.readline().rstrip()
            line = removeEmptyreads.readline().rstrip()
            rememptyprocessed = int(line.split()[0])
            line = removeEmptyreads.readline().rstrip()
            pairs = int(line.split()[0])
            line = removeEmptyreads.readline().rstrip()
            singlets = int(line.split()[0])
            
            #assert rememptyprocessed == r1ReadsCutadapt and rememptyprocessed == r2ReadsCutadapt
            #print 'this gave',pairs,'pairs that were used for mapping'
        except  (IOError, IndexError):
            rememptyprocessed = 'NA'
            pairs = 'NA'
            singlets = 'NA'
        
        try:
            bowtielines = open(path+sample+'/stderr.bowtie2.txt','r').read().split('\n')
            overallalignmentrate = float(bowtielines[-2].split()[0][:-1])/100
            totalReadsPreFilter = int(bowtielines[0].split(' ')[0])
            #assert totalReadsPreFilter == pairs
            #print overallalignmentrate,'% mapped to the genome'
        except  (IOError, IndexError):
            totalReadsPreFilter = 'NA'
            overallalignmentrate = 'NA'
        
        try:
            postmapnmarklines = open(path+sample+'/stderr.postMapNmark.txt','r').read().split('\n')
            mappedreadspassingfilters = int(postmapnmarklines[-15].split()[0])
            properlypaired = int(postmapnmarklines[-9].split()[0])/2
            #print properlypaired,' was PP',percentage(properlypaired,pairs),'%'
        except  (IOError, IndexError):
            properlypaired = 'NA'

        try:
            markDups = open(path+sample+'/'+sample+'.MarkDupsMetrix','r')
            dupsLine = markDups.read().split('\n')[7].rstrip()
            percentageDuplication   = float(dupsLine.split('\t')[7].replace(',','.'))
            markDups.close()
        except  (IOError, IndexError):
            #print sample, 'SOMETHING IS WRONG with markdups'
            percentageDuplication ='NA'
        
        #print percentageDuplication,'% duplication'
        
        try:
            hsMetrix = open(path+sample+'/'+sample+'.recalibrated.final.bam.hs_metrics.summary.txt','r')
 #       except IOError:pass
#        try:
            hsLine = hsMetrix.read().split('\n')[7].rstrip()
            onTargetPercentage      = float(hsLine.split('\t')[17].replace(',','.'))
            meanTargetCoverage      = float(hsLine.split('\t')[21].replace(',','.'))
            totalReads              = int(hsLine.split('\t')[5])
            hsMetrix.close
        except  (IOError, IndexError, ValueError):
            #print sample, 'SOMETHING IS WRONG with hsmetrics'
            onTargetPercentage = 'NA'
            meanTargetCoverage = 'NA'
            totalReads = 'NA'
        
        try:mappingRate             = float(totalReads) / float(totalReadsPreFilter*2) #float(hsLine.split('\t')[11].replace(',','.'))
        except ValueError:mappingRate = 'NA'
        #print mappingRate,'mapping rate'
        
        values = [sample, r1Reads, r1Bases+r2Bases, r1TrimmedBases+r2TrimmedBases, percentage(r1TrimmedBases+r2TrimmedBases,r1Bases+r2Bases),r1BasesCutadapt+r2BasesCutadapt, r1TrimmedBasesCutadapt+r2TrimmedBasesCutadapt,percentage(r1TrimmedBasesCutadapt+r2TrimmedBasesCutadapt,r1BasesCutadapt+r2BasesCutadapt),r1qualbases+r2qualbases,r1qualtrimbases+r2qualtrimbases,percentage(r1qualtrimbases+r2qualtrimbases,r1qualbases+r2qualbases),pairs,overallalignmentrate,mappedreadspassingfilters,percentage(mappedreadspassingfilters,pairs+pairs),properlypaired,percentage(properlypaired,pairs),percentageDuplication,onTargetPercentage,meanTargetCoverage]
        print '\t'.join([str(plupp).replace('.',',') for plupp in values])
        
        #
        # save data to dictionary
        #
        data[sample]['totalReadsPreFilter'] = totalReadsPreFilter
        data[sample]['percentageDuplication'] = percentageDuplication
        data[sample]['totalReads'] = totalReads
        data[sample]['mappingRate'] = mappingRate
        data[sample]['onTargetPercentage'] = onTargetPercentage
        data[sample]['meanTargetCoverage'] = meanTargetCoverage
        
        #
        # write to outfile
        #
        outfile.write( '\t'.join([str(plupp).replace('.',',') for plupp in values]) +'\n')
        #outfile.write( str(sample)+'\t'+str(totalReadsPreFilter*2).replace('.',',')+'\t'+str(totalReads).replace('.',',')+'\t'+str(mappingRate).replace('.',',')+'\t'+str(onTargetPercentage).replace('.',',')+'\t'+str(meanTargetCoverage).replace('.',',')+'\t'+str(percentageDuplication).replace('.',',')+'\n')
    
    outfile.close() # close outfile
    
    #
    # Make graphics
    #

def bufcount(filename):
	""" returns the number of lines in a file
	"""
	import gzip
	if filename.split('.')[-1] in ['gz','gzip']: f = gzip.open(filename)
	else: f = open(filename)
	lines = 0
	buf_size = 1024 * 1024
	read_f = f.read # loop optimization
	
	buf = read_f(buf_size)
	while buf:
            lines += buf.count('\n')
            buf = read_f(buf_size)
            f.close
	return lines

def percentage(count,total):
    if 'NA' in [total,count] or 'NANA' in [total,count]: return 'NA'
    return round(float(count) / float(total),4)
    #return round(100* float(count) / float(total),2)


if __name__ == "__main__": main()

