print ("Start")

#TODO: Watch for things with no LOC ID, or LOC IDs that look weird, including
#<Hit_def>Diplodia corticola transcription factor tfiia complex subunit mRNA</Hit_def>
#or
#<Hit_def>PREDICTED: Tursiops truncatus angiomotin like 1 (AMOTL1), transcript variant X3, mRNA</Hit_def>

import gzip, os, re, time, itertools, argparse
from sets import Set
from Bio.Blast import NCBIXML


parser = argparse.ArgumentParser(description="")
parser.add_argument('--species',nargs=1, required=True)		#length
parser.add_argument('--lengthFilter',nargs=1, required=True)		#input files
parser.add_argument('--processNames',action='store_true', required=True)
parser.add_argument('--forbiddenStrings',nargs=1)		#input files
parser.add_argument('--rmatsMode',action='store_true')
args = parser.parse_args()
species = args.species[0]
lengthFilter = int(args.lengthFilter[0])
processNames = args.processNames
forbiddenStrings = args.forbiddenStrings
rmatsMode = args.rmatsMode



class AlignmentEntry:	#Data unique to one alignment
	def __init__(self,alignment, locHit, fileNumber, extraColumns):
		self.alignment=alignment
		self.locHit=locHit
		self.fileNumber=fileNumber
		self.taxID = 0
		self.speciesName = "."
		self.extraColumns = extraColumns


class FullEntry:	#Data unique to one output stream (species, nonspecies, etc) of one file
	def __init__(self,alignmentEntries, locDict, thisLOC):
		self.alignmentEntries=alignmentEntries
		self.locDict=locDict
		self.thisLOC=thisLOC


def fileOutput(fullEntries, fileHandle,rmatsMode,extraColumnIDs):
	if(rmatsMode):
		fileHandle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("Searched ncRNA LOC ID","LOC ID of gene hit","e value","Length of sequence match","Hits on other variants","TaxID","species name","Description of gene hit","File Number"))
		for extraColumnID in extraColumnIDs:
			fileHandle.write('\t' + extraColumnID.strip('\n'))
		fileHandle.write('\n')
	else:
		fileHandle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Searched ncRNA LOC ID","LOC ID of gene hit","e value","Length of sequence match","Hits on other variants","TaxID","species name","Description of gene hit","File Number"))
	for entry in fullEntries:	#All the searches for this stream
		for alignmentEntry in entry.alignmentEntries:	#All the alignments in one search, titles live at this layer
			for hsp in alignmentEntry.alignment.hsps:	#All the high scoring matches in order of score
				displayLOC="."
				findLOC = re.search("LOC\d+",alignmentEntry.locHit)
				if(findLOC != None):	#If the hit LOC is just the title, no need to display it twice
					displayLOC = alignmentEntry.locHit
				if(rmatsMode):
					fileHandle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (entry.thisLOC,displayLOC,str(hsp.expect),str(hsp.align_length),entry.locDict[alignmentEntry.locHit]-1,str(alignmentEntry.taxID),str(alignmentEntry.speciesName),alignmentEntry.alignment.title,alignmentEntry.fileNumber))
					for extraColumn in alignmentEntry.extraColumns:
						fileHandle.write('\t' + extraColumn)
				else:
					fileHandle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entry.thisLOC,displayLOC,str(hsp.expect),str(hsp.align_length),entry.locDict[alignmentEntry.locHit]-1,str(alignmentEntry.taxID),str(alignmentEntry.speciesName),alignmentEntry.alignment.title,alignmentEntry.fileNumber))
				break	#Break to avoid printing out multiple hits per gene-to-gene match
	return

def insertTaxID(fullEntries, accDict, speciesNames, firstTime):	#Given a dictionary of accession numbers and TaxIDs, and TaxIDs and species names, attach the TaxID and species names
	for entry in fullEntries:
		for alignmentEntry in entry.alignmentEntries:
			titleParts = alignmentEntry.alignment.title.split("|")
			if not len(titleParts) == 5 and firstTime:
				print ("Strange alignment title in file %s:" % (alignmentEntry.fileNumber))
				print alignmentEntry.alignment.title
			if len(titleParts) >=3:
				if titleParts[3] in accDict:
					alignmentEntry.taxID = accDict[titleParts[3]]
					if alignmentEntry.taxID in speciesNames:
						alignmentEntry.speciesName = speciesNames[alignmentEntry.taxID]
				
	return

def finalProcessing(fullEntries, filterFileHandle):	#One more pass through, now that all the taxIDs are in place.  filterFileHandle is all the filtered results just for sanity checking
#	speciesHits
	hits = {}
	for entry in fullEntries:
		thisEntryCount = 0
		for alignmentEntry in entry.alignmentEntries:
			thisEntryCount += 1
		if(thisEntryCount != 0):
			hits[entry.thisLOC] = thisEntryCount
	return
		
#Process the accession2taxid file, append the taxID to the results
def processAcctoTax (results, speciesNames, fileHandle):
	start = True
	limit = 10000000	#Set to avoid memory issues.  Set lower if on an exceptionally low-memory machine.
	iteration = 0
	baselineTime = time.time()
	lastTime = time.time()
	accCount = 0
	accToTax = {}
	for line in fileHandle:
		if start:
			start = False
			continue
		accCount += 1
		if (time.time() - lastTime > 3):
			lastTime = time.time()
			print len(accToTax) + limit * iteration	#Just to keep track of how things are moving.  If it stalls out, set "limit" lower.
		fields = line.split("\t")
		if(len(fields) == 4):
			accToTax[fields[1]] = fields[2]
		else:
			print "Strange field in accession2taxid file:"
			print line
		if(accCount == limit * (iteration+1)):
			firstTime=False
			if iteration == 0:
				firstTime = True
			iteration += 1
			for result in results:
				insertTaxID(result, accToTax, speciesNames, firstTime)
			accToTax = {}
	insertTaxID(result, accToTax, speciesNames, firstTime)
	return

fileAppend = ""
if(lengthFilter != 0):
	fileAppend = str(lengthFilter)
with open("processThese","r") as processThese:
	firstLine = processThese.readline()
	if (len(firstLine.split('\t')) > 1):
		rmatsMode = True	#We've got extra columns, which means a few deviations in the flow, like ignoring the first line, which contains column identifiers instead of real data
		extraColumnIDs = firstLine.split('\t')[1:]
		print len(extraColumnIDs)
	else:
		processThese.seek(0)
		extraColumnIDs = []
	with open(species + fileAppend + "ncRNA","w") as speciesNCRNA:
		with open(species + fileAppend + "else","w") as speciesElse:
			with open("Non" + species + fileAppend,"w") as nonSpecies:
				speciesNCRNAResults = []
				speciesElseResults = []
				nonSpeciesResults = []
				noLocCount = 0
				for GISearched in processThese:
					extraColumns=GISearched.split('\t')[1:]	#Extra data to stick on from rMATS processing
					GISearched=GISearched.split('\t')[0].strip()
					if(os.path.isfile(os.path.dirname(os.path.realpath(__file__)) + "/data/" + GISearched.strip() + ".xml.gz")):
#						reMatch = re.match("^(\d+)\.xml\.gz$",a)
#						if not reMatch == None:
#							with gzip.open(reMatch.group(1) + ".xml.gz", "rb") as openIt:
						with gzip.open("data/" + GISearched.strip() + ".xml.gz", "rb") as openIt:
							blast_record = NCBIXML.read(openIt)
							if(rmatsMode):
								LOCtouse = extraColumns[0]
							else:
								reLOC = re.search("LOC\d+",blast_record.query)
								if(reLOC == None):
									print "Error, LOC not found in blast record query:"
									print blast_record.query
									continue
								LOCtouse = reLOC.group(0)

							speciesNCRNAAlignments = []
							speciesElseAlignments = []
							nonSpeciesAlignments = []
							speciesNCRNALocDict = {}
							speciesElseLocDict = {}
							nonSpeciesLocDict = {}
							for alignment in blast_record.alignments:
								thisLocHit = re.search("LOC\d+",alignment.title)
								if(alignment.hsps[0].align_length < lengthFilter):
									continue
								if(thisLocHit == None):	#No included LOC, so match the whole title field.  This isn't perfect, as sometimes variants aren't labeled identically.
									thisLocHit = alignment.title
									noLocCount += 1
								else:
									thisLocHit = thisLocHit.group(0)
								forbidden = False
								for forbiddenString in forbiddenStrings:
									if forbiddenString in alignment.title:
										forbidden = True
										break
								if(forbidden):
									continue
								if(species in alignment.hit_def):
									if("ncRNA" in alignment.hit_def and (not str(LOCtouse) in alignment.hit_def)):		#ncRNA
										if(not (thisLocHit in speciesNCRNALocDict.keys())):
											speciesNCRNALocDict[thisLocHit] = 1
											speciesNCRNAAlignments.append(AlignmentEntry(alignment,thisLocHit,GISearched,extraColumns))
										else:
											speciesNCRNALocDict[thisLocHit] += 1
									if((not "ncRNA" in alignment.hit_def) and (not str(LOCtouse) in alignment.hit_def)):	#not ncRNA
										if(not (thisLocHit in speciesElseLocDict.keys())):
											speciesElseLocDict[thisLocHit] = 1
											speciesElseAlignments.append(AlignmentEntry(alignment,thisLocHit,GISearched,extraColumns))
										else:
											speciesElseLocDict[thisLocHit] += 1
								else:		#Non-species
									if(not (thisLocHit in nonSpeciesLocDict.keys())):
										nonSpeciesLocDict[thisLocHit] = 1
										nonSpeciesAlignments.append(AlignmentEntry(alignment,thisLocHit,GISearched,extraColumns))
									else:
										nonSpeciesLocDict[thisLocHit] += 1
							speciesNCRNAResults.append(FullEntry(speciesNCRNAAlignments,speciesNCRNALocDict,LOCtouse))
							speciesElseResults.append(FullEntry(speciesElseAlignments,speciesElseLocDict,LOCtouse))
							nonSpeciesResults.append(FullEntry(nonSpeciesAlignments,nonSpeciesLocDict,LOCtouse))
					else:
						print "Missing BLAST result file: " + GISearched + ".xml.gz"
					

				print "No LOCs: " + str(noLocCount)
				speciesNCRNAResults.sort(key=lambda results: len(results.alignmentEntries), reverse=True)
				nonSpeciesResults.sort(key=lambda results: len(results.alignmentEntries), reverse=True)
				speciesElseResults.sort(key=lambda results: len(results.alignmentEntries), reverse=True)
				if(processNames):
					accToTax = {}
					start = True
					limit = 10000000
					accCount = 0
					iteration = 0
					speciesNames = {}
					with open("names.dmp","r") as speciesNamesIn:
						for line in speciesNamesIn:
							#nameFields = line.split("|")
							nameFields = [field.strip("\t") for field in line.split("|")]
							if "scientific name" in nameFields:
								if len(nameFields) >= 2:
									speciesNames[nameFields[0]] = nameFields[1]
								else:
									print "Strangely formatted line in names.dmp"
									print line
					print "reading gb accession table"
					with open("nucl_gb.accession2taxid","r") as accToTaxFile:
						processAcctoTax([speciesNCRNAResults,speciesElseResults,nonSpeciesResults],speciesNames,accToTaxFile)

				with open("filteredResults","w") as filteredResults:
					finalProcessing(speciesNCRNAResults, filteredResults)
					finalProcessing(nonSpeciesResults, filteredResults)
					finalProcessing(speciesElseResults, filteredResults)
				fileOutput(speciesNCRNAResults,speciesNCRNA,rmatsMode,extraColumnIDs)
				fileOutput(nonSpeciesResults,nonSpecies,rmatsMode,extraColumnIDs)
				fileOutput(speciesElseResults,speciesElse,rmatsMode,extraColumnIDs)
				nonSpecies.write("\n")
				speciesElse.write("\n")
				speciesNCRNA.write("\n")
		
	

print ("End")


