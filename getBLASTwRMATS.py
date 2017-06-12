print ("Start")

import re
import time
import gzip
from Bio.Blast import NCBIWWW
import os.path
import sys
from pyfasta import Fasta

class NCRNAEntry:
	def __init__(self,LOCID,coord,GI,length,numExons):
		self.LOCID=LOCID
		self.coord=coord
		self.GI=GI
		self.length=length
		self.numExons=numExons

class RMATSEntry:
	def __init__(self,geneSymbol, chromosome, startCoord, endCoord, asType, rmatsID, pValue, FDR, mxeB = False):
		self.geneSymbol=geneSymbol
		self.chromosome = chromosome
		self.startCoord = startCoord
		self.endCoord = endCoord
		self.asType = asType
		self.rmatsID = rmatsID
		self.startCoord = startCoord
		self.endCoord = endCoord
		self.pValue = pValue
		self.FDR = FDR
		self.mxeB = mxeB	#Is it the second part of an MXE?

def getBLASTseq(seq,filename):
        if(os.path.isfile(filename)):
                return
	print ("Getting BLAST result" + filename)
	startTime = time.time()
	result_handle = NCBIWWW.qblast("blastn","refseq_rna",seq,hitlist_size=500,megablast=True)	#All RNA now
	print ("Got BLAST result in %d seconds" % (time.time() - startTime))
	with gzip.open(filename, "wb") as saveIt:
		saveIt.write(result_handle.read())
	result_handle.close()
	return

#def getBLASTseq2(seq,filename):
#	print seq + " " + filename
#	return

withGIs = True
withCoords = True
filterPValue = 0.005
filterFDR = 0.05

if len(sys.argv) < 2:
	print "Error, give path to rMATS results dir"
	exit()
dirToRMatsResults = sys.argv[1]
rmatsResultFiles = {	"A3SS" : dirToRMatsResults + "/MATS_output/" + "A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt",
			"A5SS" : dirToRMatsResults + "/MATS_output/" + "A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt",
			"MXE" : dirToRMatsResults + "/MATS_output/" + "MXE.MATS.ReadsOnTargetAndJunctionCounts.txt",
			"RI" : dirToRMatsResults + "/MATS_output/" + "RI.MATS.ReadsOnTargetAndJunctionCounts.txt",
			"SE" : dirToRMatsResults + "/MATS_output/" + "SE.MATS.ReadsOnTargetAndJunctionCounts.txt" }
for ASType in rmatsResultFiles:
	if(not os.path.isfile(rmatsResultFiles[ASType])):
		print "rMATS results not found :" + rmatsResultFiles[ASType]
		exit()
rmatsResults = []

for ASType in rmatsResultFiles:
	pValueField = 0
	FDRField = 0
	with open(rmatsResultFiles[ASType], "r") as resultIn:
		firstLine = resultIn.readline().split("\t")
		for i in xrange(len(firstLine)):
			if firstLine[i] == "PValue":
				pValueField = i
			if firstLine[i] == "FDR":
				FDRfield = i
		if(pValueField == 0 or FDRfield == 0): 
			print "Constraints not found, check for change in rMATS output format"
		for line in resultIn:
			fields = line.split("\t")
			if len(fields) < 7:
				print "Strangely formatted input, not enough fields" + line
				continue
			if(float(fields[pValueField]) >= filterPValue or float(fields[FDRfield]) >= filterFDR):
				continue
			rmatsResults.append(RMATSEntry(fields[2].strip("\""),fields[3],int(fields[5]),int(fields[6]),ASType,fields[0],fields[pValueField],fields[FDRfield]))
			if(ASType == "MXE"):	
				rmatsResults.append(RMATSEntry(fields[2].strip("\""),fields[3],int(fields[7]),int(fields[8]),ASType,fields[0],fields[pValueField],fields[FDRfield], True))	#Get the second exon for MXE too
print len(rmatsResults)



#TODO: This mode thing is getting dumb, probably take it out and just do everything, outputting only as necessary
with open("GCF_001433935.1_IRGSP-1.0_genomic.gbff", "r") as gbffIn:
        
	output, numExons, thisGI, thisLOCID, foundOne = [], 0, "", "", False
	entries = []
	mode = 0	#Mode: 0 to findNC, 2 then 1 to find gene and GI.  Will behave badly if either are ever absent, but will work okay if they're in either order
	for line in gbffIn:
		if mode == 0:
			if "ncRNA" in line and not("Features" in line or "/ncRNA" in line):
				matchThis = line.strip()
				while len(re.findall('\(',matchThis)) != len(re.findall('\)',matchThis)):
					matchThis += next(gbffIn).strip()
				exonCoords = re.finditer('(\d+)\.\.(\d+)',matchThis)	#This discards e.g. < signs, so "<150..300" is treated as "150..300"
				length = 0
				numExons = 0
				for eCoord in exonCoords:
					numExons += 1
					length += int(eCoord.group(2)) - int(eCoord.group(1))
				mode = 2
				foundOne = True

		if mode != 0:
			if "/gene" in line:
				result = re.search('LOC[0-9]+',line)
				if result:
					thisLOCID = result.group(0)
				else:
					print "Error, bad LOC format"
				if(not withGIs):
					mode = 0
				else:
					mode -= 1
			if ("/db_xref=\"GI:" in line) and withGIs:
				result = re.search('GI:[0-9]+',line)
				if result:
					thisGI=result.group(0).lstrip("GI:")
				else:
					print "Error, bad GI format"
				mode -= 1
		if mode == 0 and foundOne:	#About to continue

			#	def __init__(LOCID,coord,GI,length,numExons):
			entries.append(NCRNAEntry(thisLOCID,matchThis.lstrip("ncRNA").strip(),thisGI,length,numExons))
			result = "%s with length %d and %d exons" % (thisLOCID,length,numExons)
			if(withGIs):
				result = "%s GI: %s" % (result,thisGI)
			if(withCoords):
				result = "%s : %s" % (result,matchThis.lstrip("ncRNA").strip())
			output.append(result)
			foundOne = False


#	def __init__(self,geneSymbol, chromosome, startCoord, endCoord, asType, rmatsID, pValue, FDR, mxeB = False)
#def getBLASTseq(seq,filename)
count = 0
#genome = Fasta('/media/alex/Two/BioinfoFiles/Oryza_sativa_japonica_Ensembl_IRGSP-1.0/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/WholeGenomeFasta/genome.fa')
genome = Fasta('y:/BioinfoFiles/Oryza_sativa_japonica_Ensembl_IRGSP-1.0/Oryza_sativa_japonica/Ensembl/IRGSP-1.0/Sequence/WholeGenomeFasta/genome.fa')
#print genome.keys()
#print genome['5'][601370:601402]
with open("processThese",'w') as processThese:
	processThese.write(".\tLOC ID\tPValue\tFDR\tchr\tstartCoord\tendCoord\tASType\n")
	with open('matches','w') as matchesOut:
		matchesOut.write("Gene symbol\tAS type\tID\n")
		for rmatsResult in rmatsResults:
			for entry in entries:
				if rmatsResult.geneSymbol == entry.LOCID:
					count += 1
					matchesOut.write(entry.LOCID + "\t" + rmatsResult.asType + "\t" + rmatsResult.rmatsID + "\n")
					if(rmatsResult.mxeB):
						getBLASTseq(genome[rmatsResult.chromosome][int(rmatsResult.startCoord):int(rmatsResult.endCoord)],rmatsResult.asType + rmatsResult.rmatsID + "B.xml.gz")
						processThese.write(rmatsResult.asType + rmatsResult.rmatsID + "B\t" + rmatsResult.geneSymbol + '\t' + rmatsResult.pValue + '\t' + rmatsResult.FDR + '\t' + rmatsResult.chromosome + '\t' + str(rmatsResult.startCoord) + '\t' + str(rmatsResult.endCoord) + '\t' + rmatsResult.asType + '\n')
					else:
						getBLASTseq(genome[rmatsResult.chromosome][int(rmatsResult.startCoord):int(rmatsResult.endCoord)],rmatsResult.asType + rmatsResult.rmatsID + ".xml.gz")
						processThese.write(rmatsResult.asType + rmatsResult.rmatsID + '\t' + rmatsResult.geneSymbol + '\t' + rmatsResult.pValue + '\t' + rmatsResult.FDR + '\t' + rmatsResult.chromosome + '\t' + str(rmatsResult.startCoord) + '\t' + str(rmatsResult.endCoord) + '\t' + rmatsResult.asType + '\n')
					break

print count
print len(output)




