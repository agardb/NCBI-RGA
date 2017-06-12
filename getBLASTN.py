print ("Start")

import re
import time
import gzip
from Bio.Blast import NCBIWWW
import os.path

class Entry:
	def __init__(self,LOCID,coord,GI,length,numExons):
		self.LOCID=LOCID
		self.coord=coord
		self.GI=GI
		self.length=length
		self.numExons=numExons

def getBLAST(GI):	#
	print ("Getting BLAST result")
	startTime = time.time()
	result_handle = NCBIWWW.qblast("blastn","refseq_rna",GI,hitlist_size=500,megablast=True)	#All RNA now
#	result_handle = NCBIWWW.qblast("blastn","refseq_rna",GItouse,db_genetic_code=39947,hitlist_size=5,megablast=True)
	print ("Got BLAST result in %d seconds" % (time.time() - startTime))
	with gzip.open(str(GI) + ".xml.gz", "wb") as saveIt:
		saveIt.write(result_handle.read())
	result_handle.close()
	return


withGIs = True
withCoords = False
justCount = False
onlyLongest = True
onlyShortest = False
if onlyLongest or onlyShortest:
	withCoords = True
count = 0
#TODO: This mode thing is getting dumb, probably take it out and just do everything, outputting only as necessary
#with open("/media/alex/Two/BioinfoFiles/GCF_001433935.1_IRGSP-1.0_genomic.gbff", "r") as gbffIn:
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
				if justCount:
					count += 1
					mode = 0

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
			entries.append(Entry(thisLOCID,matchThis.lstrip("ncRNA").strip(),thisGI,length,numExons))
			result = "%s with length %d and %d exons" % (thisLOCID,length,numExons)
			if(withGIs):
				result = "%s GI: %s" % (result,thisGI)
			if(withCoords):
				result = "%s : %s" % (result,matchThis.lstrip("ncRNA").strip())
			output.append(result)
			foundOne = False


if(justCount):
	print count
output.sort()

with open('ncRNAout', 'w') as outFile :
	outFile.truncate()
	for a in output:
		print>>outFile,a



print len(output)





entries.sort(key=lambda entry: entry.LOCID)
#def __init__(self,LOCID,coord,GI,length,numExons):
with open('processThese','w') as processThese:
	if onlyLongest:
		if(len(entries) > 0):
		        longestFound = entries[0]
		for i in xrange(0,len(entries)):
		        #print str(entries[i].LOCID) + " " + str(entries[i].GI) + " " + str(entries[i].length) + " Best: " + str(longestFound.GI)
		        if(entries[i].LOCID == longestFound.LOCID):
		                if(entries[i].length > longestFound.length):
		                        longestFound = entries[i]
		        else:		
		                if(not os.path.isfile(os.path.dirname(os.path.realpath(__file__)) + "/" + longestFound.GI + ".xml.gz")):
		                        getBLAST(longestFound.GI)
				processThese.write(longestFound.GI + '\n')
		                longestFound = entries[i]
		if(longestFound.GI != None):
		        if(not os.path.isfile(os.path.dirname(os.path.realpath(__file__)) + "/" + longestFound.GI + ".xml.gz")):
		                getBLAST(longestFound.GI)
			#processThese.write(longestFound.GI + '\n')

	if onlyShortest:
		if(len(entries) > 0):
		        shortestFound = entries[0]
		for i in xrange(0,len(entries)):
		        #print str(entries[i].LOCID) + " " + str(entries[i].GI) + " " + str(entries[i].length) + " Best: " + str(shortestFound.GI)
		        if(entries[i].LOCID == shortestFound.LOCID):
		                if(entries[i].length < shortestFound.length):
		                        shortestFound = entries[i]
		        else:		
		                if(not os.path.isfile(os.path.dirname(os.path.realpath(__file__)) + "/" + shortestFound.GI + ".xml.gz")):
		                        getBLAST(shortestFound.GI)
				processThese.write(longestFound.GI + '\n')
		                shortestFound = entries[i]
		if(shortestFound.GI != None):
		        if(not os.path.isfile(os.path.dirname(os.path.realpath(__file__)) + "/" + shortestFound.GI + ".xml.gz")):
		                getBLAST(shortestFound.GI)
			#processThese.write(longestFound.GI + '\n')






