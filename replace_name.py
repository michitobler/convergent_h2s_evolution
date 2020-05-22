#!/usr/bin/env python
# Author: Joanna L. Kelley
# replace_name.py
# Replaces name in fasta with only one fasta entry with new name

import sys
from optparse import  OptionParser
import random
import os


USAGE = """USAGE: replace_name.py --f <fasta file> --r <text> --o <output filename> """

parser = OptionParser(USAGE)
parser.add_option("--f",dest="fastaFile",help="name of fasta file to be processed")
parser.add_option("--o",dest="outFile",help="name of fasta file to be written as outfile")
parser.add_option("--r",dest="textData",help="text to replace fasta header")

(options,args)=parser.parse_args()


if options.fastaFile is None:
    parser.error('fasta file not given')
if options.outFile is None:
    parser.error('output filename not given')
if options.textData is None:
    parser.error('text for replacement not given')
    
    
    
fasta = open(options.fastaFile,'r')
newname = options.textData
newfasta = open(options.outFile, 'w')
count = 1

for line in fasta:
    if line.startswith('>'):
    	tofill = '>' + newname + '\n'
    	newfasta.write(tofill)
    	newname = newname + str(count)
    	count += 1		
    else:
        newfasta.write(line)

newfasta.write('\n')

fasta.close()
newfasta.close()

