import argparse
import sys

#Commandline arguments instruction.
parser  = argparse.ArgumentParser(prog="get_data.py",
  description = "Our simulation outputs a file that can't be used with gnuplot. This is intended; this file is able to convert.")  
parser.add_argument('-n', '--number', help='Which particular data set is wanted?.', default = 0, action='store', type = int)   
parser.add_argument('-f', '--file', help='What data file?.', default = "rough_sweep.txt", action='store', type = str)   
args    = parser.parse_args() 


number  = args.number
filename = args.file
reader    = open(filename, 'r')

lines = reader.readlines()


dataset = -1;

for line in lines:
        line = line.rstrip()
        
        if "Report" in line:
                dataset += 1
                if dataset == number:
                        print >> sys.stderr, "Printing dataset %d" % dataset
        elif dataset == number:
                print line
        
reader.close()
