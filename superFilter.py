#Filters blast results by various input filters

#inputFiles = ["Oryza sativa JaponicancRNA100.txt","Oryza sativa Japonicaelse100.txt","NonOryza sativa Japonica100.txt","Oryza sativa JaponicancRNA100rMATS.txt","Oryza sativa Japonicaelse100rMATS.txt","NonOryza sativa Japonica100rMATS.txt"]
#filterFiles = ["blastFinalEandLenFilter","blastFinalNoFilter","blastFinalEfilter","blastFinalLenFilter"]


import sys, os.path, math, re, argparse

parser = argparse.ArgumentParser(description="Filters multiple blast results by multiple input filters")
parser.add_argument('-i',nargs='*')
parser.add_argument('-f',nargs='*')
args = parser.parse_args()
filterFiles = args.f
inputFiles = args.i


filterLOCs = []
for filterStyle in filterFiles:
	with open(filterStyle,"r") as inputFilter:
		filterLOCs = []
		for line in inputFilter:
			parseIt = re.match(r"(\S+)	.*", line)
			filterLOCs += [parseIt.group(1)]
	#print(len(filterLOCs))
	for inputFile in inputFiles:
		with open(inputFile,"r") as inputFileOpen:
			with open(filterStyle + inputFile,"w") as output:
				removedLines = 0
				top = True
				for line in inputFileOpen:
					if(top):
						output.write(line)
						top = False
					else:
						parseIt = re.match(r"(\S+)	(\S+)	.*", line)
						if(parseIt.group(1) in filterLOCs or parseIt.group(2) in filterLOCs):
							removedLines += 1
						else:
							output.write(line)
				#print removedLines

print("Done")
