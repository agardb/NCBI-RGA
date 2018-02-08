#Filters blast result files based on length

import sys, argparse

parser = argparse.ArgumentParser(description="Filters blast result files based on length")
parser.add_argument('-l',nargs='*')	#length
parser.add_argument('-i',nargs='*')	#input files
args = parser.parse_args()
lengths = args.l
inputFiles = args.i

def filterFile (fileIn, fileOut, length):
	for line in fileIn:
		fields = line.split("\t")
		try:
			if(len(fields) >= 4 and int(fields[3]) >= length):
				fileOut.write(line)
		except ValueError:
			fileOut.write(line)
	return


for length in lengths:
	for inputFile in inputFiles:
		with open(inputFile,"r") as fileIn:
			with open(inputFile+str(length),"w") as fileOut:
				filterFile(fileIn,fileOut,int(length))
