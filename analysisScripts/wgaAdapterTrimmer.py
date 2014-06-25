#!/usr/bin/env python

import sys

sys.stderr.write('Running wgaAdapterTrimmer.py\n')

ADAPTERSEQUENCE = 'TGTGTTGGGTGTGTTTGG'

def main():
    indata = getComandLineOptions()

    TOTALBASES = 0
    TRIMMEDBASES = 0
    READSPROCESSED = 0
    
    for read in readGenerator(indata.infile):
        TOTALBASES,TRIMMEDBASES,read = trimRead(read, indata.maxDist,TOTALBASES,TRIMMEDBASES)
        printRead(read)
        READSPROCESSED += 1
    
    sys.stderr.write('Processed a total of\t'+str(READSPROCESSED)+'\treads. ('+indata.infile+')\n')
    sys.stderr.write('Processed a total of\t'+str(TOTALBASES)+'\tbases ('+indata.infile+').\n')
    sys.stderr.write('trimmed a total of\t'+str(TRIMMEDBASES)+'\tbases in the start of reads ('+indata.infile+').\n')
    sys.stderr.write('wgaAdapterTrimmer.py done exiting ...\n')

def trimRead(read,maxDist,TOTALBASES,TRIMMEDBASES):
    import sys
    header, sequence, qual = read
    assert len(sequence) == len(qual), '\n\nError: sequence and qual has different lengths.\n\n'
    dist = levenshtein(sequence[:len(ADAPTERSEQUENCE)], ADAPTERSEQUENCE)
    TOTALBASES += len(sequence)
    if float(dist) <= maxDist:
        #sys.stderr.write('Trimming '+header+' edit distance '+str(dist)+' <= '+str(float(len(ADAPTERSEQUENCE))*0.2)+'.\n')
        sequence = sequence[30:]
        qual = qual[30:]
        TRIMMEDBASES += 30
    
    assert len(sequence) == len(qual), '\n\nError: sequence and qual has different lengths.\n\n'
    return TOTALBASES,TRIMMEDBASES,[header, sequence, qual]

def printRead(read):
    import sys
    header, sequence, qual = read
    assert len(sequence) == len(qual), '\n\nError: sequence and qual has different lengths.\n\n'
    sys.stdout.write(header+'\n'+sequence+'\n+\n'+qual+'\n')

def readGenerator(fastq1):
    
    #
    # imports
    #
    import gzip

    if fastq1.split('.')[-1] in ['gz','gzip']: file1 = gzip.open(fastq1)
    else: file1 = open(fastq1,'r')
     
    while 'NOT EOFError':
        try:
            header1   = file1.readline().rstrip()
            sequence1 = file1.readline().rstrip()
            trash     = file1.readline().rstrip()
            qual1     = file1.readline().rstrip()
            if not header1: break
            yield header1, sequence1, qual1
        except EOFError: break

def getComandLineOptions():
    """ function that gets the indata from the commandline """

    import argparse
    import sys
    
    argparser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description='Trimmes the first 30 bases of reads where '+ADAPTERSEQUENCE+' is found.')

    # All programs
    argparser.add_argument('-i',dest='infile',	metavar='<fastq-file>',type=str,required=True, help='Indata "fastq"-file.')
    argparser.add_argument('-e', dest='maxDist', metavar='N',type=float, required=False,default=float(len(ADAPTERSEQUENCE))*0.2, help='Max edit distance to '+ADAPTERSEQUENCE+' (levenshtein, default float(len(ADAPTERSEQUENCE))*0.2).')
    indata = argparser.parse_args(sys.argv[1:])
    
    return indata

def hamming_distance(s1, s2):
	assert len(s1) == len(s2), 'Error: '+str(len(s1)) + ' != ' + str(len(s2))
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def levenshtein(s1, s2):
	if len(s1) < len(s2):
		return levenshtein(s2, s1)
	if not s1:
		return len(s2)
	previous_row = xrange(len(s2) + 1)
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
		for j, c2 in enumerate(s2):
			insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
			deletions = current_row[j] + 1       # than s2
			substitutions = previous_row[j] + (c1 != c2)
			if c1 == 'N' or c2 == 'N': substitutions -= 1 #if N then no mismatch
			current_row.append(min(insertions, deletions, substitutions))
		previous_row = current_row
	return previous_row[-1]

if __name__ == "__main__":
    main()
