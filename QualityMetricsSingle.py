#IMPORTS
import sys
import re
import csv
import shutil
import os
import itertools
#System Argument IMPORTS
annovarID = sys.argv[1]
ALSconsortiumFilePath = sys.argv[2]
dateHolder = sys.argv[3]
#VaribleSpecificToVCFFormat
header ="#"
spacebeforeDivider= "	"
startofdata = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
#VARIABLE/Array Definitions
newdirectory = ALSconsortiumFilePath + "ALS_consortium/" + dateHolder + "_annotated/Original_Annotation"
if os.path.isdir(newdirectory)==False:
    os.mkdir(newdirectory)#
copylocation = ALSconsortiumFilePath + "ALS_consortium/" + dateHolder + "_annotated/Original_Annotation/hg38_CGND-HDA"+ "-" + annovarID + ".filtered_acmg_als.avinput_updated.hg38_multianno.csv"#
annotatedvcf = ALSconsortiumFilePath + "ALS_consortium/" + dateHolder + "_annotated/hg38_CGND-HDA"+ "-" + annovarID + ".filtered_acmg_als.avinput_updated.hg38_multianno.csv"#
annotatedvcffull = ALSconsortiumFilePath + "ALS_consortium/" + dateHolder + "_annotated/hg38_CGND-HDA"+ "-" + annovarID + ".filtered_acmg_als.avinput_updated_Metrics.hg38_multianno.csv"#
holderVCF= "CGND-HDA-"+annovarID+".filtered_acmg_als.vcf"
avInputPath = ALSconsortiumFilePath + "ALS_consortium/AnnovarReady/" + dateHolder + "_annovar_ready/CGND-HDA"+ "-" + annovarID + ".filtered_acmg_als.avinput"
l = open(avInputPath, 'r')
messageAvInput = l.read()
ArrayAvInput = messageAvInput.splitlines()
alignment = []
my_filePATH = ALSconsortiumFilePath + "ALS_consortium/" + dateHolder + "/CGND-HDA"+ "-" + annovarID + ".filtered_acmg_als.vcf"
f = open(my_filePATH,'r')
messageVCF = f.read()
ArrayVCF = messageVCF.splitlines()
spaces=[]
metricscell = ""
metricsheadercell = ""
QUAL = ""
FILTERRESULT = ""
chrholder= ""
posIDholder= ""
lineholder = ""
extracted = []
DataRichArrayVCF = []
counter = 0
valcounter = 0
datafound = 0
splits = []
vcfmatchpairs = []
adjustment = 0
holderadjustment = 0
#Class Definitions
class qualitymetrics:
        def __init__(self, QUAL, FILTER, Metrics, MetricsHeader, Chr, PosID):
            self.QUAL = QUAL
            self.FILTER = FILTER
            self.Metrics = Metrics
            self.MetricsHeader = MetricsHeader
            self.Chr = Chr
            self.PosID = PosID
extracted.append(qualitymetrics("QUAL", "FILTER", "Genotype Field", "GT Field Header", "CHR", "PosID")) #Header

for i in ArrayVCF:#separating data from header
    if datafound == 1:
        DataRichArrayVCF.append(i)
    if datafound == 0:
        if (i.find(startofdata) != -1):
            datafound = 1
for lineholder in DataRichArrayVCF:#extract all data from each line i and place in array called extracted using the class qualitymetrics
    spaces=[m.start() for m in re.finditer(spacebeforeDivider, lineholder)]
    chrholder = lineholder[0:(spaces[0])]
    posIDholder = lineholder[(spaces[0]+1):(spaces[1])]
    QUAL = lineholder[(spaces[4]+1):(spaces[5])]
    FILTERRESULT = lineholder[(spaces[5]+1):(spaces[6])]
    metricsheadercell = lineholder[(spaces[7]+1):(spaces[8])]
    metricscell = lineholder[(spaces[8]+1):]
    extracted.append(qualitymetrics(QUAL, FILTERRESULT, metricscell, metricsheadercell, chrholder, posIDholder))
shutil.move(annotatedvcf, copylocation)
for j in ArrayAvInput:#extracting all data from avinput file inorder to allign the quality metric values appropriatly
    spaces=[o.start() for o in re.finditer(spacebeforeDivider, j)]
    chrholder = j[0:(spaces[0])]
    posIDholder = j[(spaces[0]+1):(spaces[2])]
    QUAL = j[(spaces[5]+1):(spaces[6])]
    FILTERRESULT = "N/A"
    metricscell = "N/A"
    metricsheadercell = "N/A"
    alignment.append(qualitymetrics(QUAL, FILTERRESULT, metricscell, metricsheadercell, chrholder, posIDholder))
for s in extracted:#find the values in extracted that have matched pairs that have the same QUAL value and are adjacent to each other
    if s.QUAL == extracted[valcounter+1].QUAL:
        vcfmatchpairs.append(valcounter)
    valcounter = valcounter +1
    if valcounter == (len(extracted)-1):
        break
datafound = 0
if len(alignment) != (len(extracted)-1):#Appends/Inserts repeats of the objects in the extracted list to make sure that the values are properly aligned
    counter = 0
    for n in alignment:
        if len(vcfmatchpairs) != 0:
            for x in vcfmatchpairs:
                holderadjustment = int(counter) + 1 - int(adjustment)
                if (int(x) == holderadjustment):#calibrate this number
                    datafound = 1
        if ((n.QUAL == alignment[len(alignment)-2].QUAL) and (datafound == 0) and (n.QUAL == alignment[counter+1].QUAL)):
            extracted.append(qualitymetrics(extracted[counter+1].QUAL, extracted[counter+1].FILTER, extracted[counter+1].Metrics, extracted[counter+1].MetricsHeader, extracted[counter+1].Chr, extracted[counter+1].PosID))
            datafound = 1
            adjustment = adjustment + 1
        if ((n.QUAL == alignment[counter+1].QUAL) and (datafound == 0)):
            extracted.insert(counter+2, qualitymetrics(extracted[counter+1].QUAL, extracted[counter+1].FILTER, extracted[counter+1].Metrics, extracted[counter+1].MetricsHeader, extracted[counter+1].Chr, extracted[counter+1].PosID))
            adjustment = adjustment + 1
        counter = counter + 1
        datafound = 0
        if counter == (len(alignment)-1):
            break
    if len(alignment) != (len(extracted)-1):
        print("ERROR: ALIGNMENT UNSUCCESSFUL")
        datafound = 1
#csv backup and adjustment
with open(copylocation, "rb") as in_file:#if this needs to be in binary convert the r to rb  remove   , newline=''
    csv_file_in = csv.reader(in_file, dialect='excel')
    with open(annotatedvcffull,"wb") as out_file:#if this needs to be in binary convert the w to wb  remove   , newline=''
        csv_file_out = csv.writer(out_file, dialect='excel')
        counter = 0
        for row in csv_file_in:
            if datafound == 1:
                break
            if datafound == 0:
                csv_file_out.writerow(row[0:3]+ [extracted[counter].QUAL, extracted[counter].FILTER, extracted[counter].MetricsHeader, extracted[counter].Metrics] + row[3:]) #Add obj call for chr and pos id if needed for test
            counter = counter + 1
if datafound == 0:
    print("QualityMetricSingle.py Run Succesfully With "+ annovarID +" ID")
