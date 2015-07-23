# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 15:37:14 2015

@author: Jason Wang
@title: qPCR Automated Analysis
"""
import math 

def sampleNameAnalysis(list):
    sampleName = list[1]
    indexUnknown = list.index("Unknown")
    originalWhiteSpaces = 0
    if indexUnknown > 4:
        whiteSpaces = indexUnknown - 4
        originalWhiteSpaces = whiteSpaces
        sampleName = list[1]
        next = 1
        while whiteSpaces > 0: 
            sampleName = sampleName + "_" + str(list[1+next])
            whiteSpaces -= 1
            next += 1
    return sampleName, originalWhiteSpaces

count = 0
F1 = open("20150513_qPCRYChrom.txt", "rU") #direct output from the qPCR machine

#skip first 11 header lines
while count < 11:
    F1.readline()
    count += 1

sampleList = []
medianCtList = []
primerList = []

for line in F1:
    list = line.strip().split() 
    if len(list) > 0 and list[0].isdigit() == True: #only sample lines are analyzed
        #Obtaining Sample Name...        
        sampleName, whiteSpaces = sampleNameAnalysis(list)
        print sampleName, whiteSpaces
        sampleList.append(sampleName)
        
        #Obtaining Median Cycle Threshold...
        #the median Ct comes right after the first instance of Tm"
        tmIndex = list.index("Tm")
        medianCt = list[tmIndex + 1]
        if medianCt.find(".") == -1: #no decimal point found
            medianCtList.append("NA")
        else:
            medianCtList.append(medianCt)
        
        #Obtaining amplicon name
        #Using "Unknown" as an index marker from the "Task" column
        unknownIndex = list.index("Unknown")
        detectorEndIndex = unknownIndex - 1
        detectorName = ""
        for index in xrange(2 + whiteSpaces, detectorEndIndex):
            detectorName = detectorName + list[index]
        primerList.append(detectorName)

F1.close()

uniquePrimers = []
for primer in primerList:
    if primer not in uniquePrimers:
        uniquePrimers.append(primer)
        
print uniquePrimers

#Spot Check:      
#print medianCtList, len(medianCtList)
#print sampleList, len(sampleList)
#print primerList, len(primerList)
      
#Create a dictionary with the key: sample name and value: [median Ct, primer]
masterSampleDictionary = {}

#assuming medianCtList and sampleNameList have the same length
for index in xrange(0, len(medianCtList), 3): #skip by 3 since experiments are run in triplicate
    masterSampleDictionary[primerList[index] + " " + sampleList[index]] = [medianCtList[index], primerList[index]]

print masterSampleDictionary

def primerAnalysis(masterSampleDictionary, startIndex, F2):
    print "START INDEX: " + str(startIndex)
    sortedSample = sorted(masterSampleDictionary)
    primer = sortedSample[startIndex].strip().split(" ")[0]
    print "Current Primer Undergoing Analysis: " + primer
    numSamples = 0
    while True:   
            try:
                numSamples = int(raw_input("How many samples were evaluated using primer " + primer + "?: "))
                break #breaks out of while loop if no error is encountered
            except ValueError:
                print "The value type you entered was invalid. Try again!"
    
    numStandards = 0
    while True:   
        try:
            numStandards = int(raw_input("How many standards were evaluated using primer " + primer + "?: (e.g. 0, 5, etc.) "))
            break #breaks out of while loop if no error is encountered
        except ValueError: 
            print "The value type you entered was invalid. Try again!"

    dilutionFactor = -1
    if numStandards > 0:
        while True:
            try:
                dilutionFactor = raw_input("By what factor did you dilute the standards?: ")
                break #breaks out of while loop if no error is encountered
            except ValueError: 
                print "The value type you entered was invalid. Try again!"
        
    #RE-PROMPT TO MAKE SURE USER INPUTS CORRECT VALUE AND AN INTEGER VALUE
    totalSamples = numSamples + numStandards
    
    #TEST: ALK has 0 standards and 8 samples
    sampleDictionary = {}
    standardDictionary = {}
    for index in xrange(startIndex, startIndex + totalSamples): 
        print sortedSample[index]
        print index
        if "Standard" in sortedSample[index] or "standard" in sortedSample[index]:
            standardDictionary[sortedSample[index]] = masterSampleDictionary[sortedSample[index]]
        else:
            sampleDictionary[sortedSample[index]] = masterSampleDictionary[sortedSample[index]]
    
    print sampleDictionary
    print standardDictionary
    
    #calculate standard Ct Differences
    standardCtDifferences = {}
    
    if numStandards == 0: #no standards
        standardCtDifferences = "NA"
    else: 
        standardSorted = sorted(standardDictionary)
        for num in xrange(1, int(numStandards)):
            key = "S" + str(num + 1) + " - S" + str(num)
            standard1 = standardSorted[num - 1]
            standard2 = standardSorted[num]
            standard1Ct = float(standardDictionary[standard1][0])
            standard2Ct = float(standardDictionary[standard2][0])
            value = round(standard2Ct - standard1Ct, 6)
            standardCtDifferences[key] = value
    
    print "Standard Ct Difference: " 
    if numStandards == 0:
        print standardCtDifferences
    else:
        for key in standardCtDifferences:
            print key + ": " + str(standardCtDifferences[key]) 
    
    #calculate primer efficiencies
    efficiency = 0.0 #temporary initalization
    if numStandards == 0:
        print "Your qPCR reaction did not include standards for this primer."
        while True:   
            try:
                efficiency = float(raw_input("Please input the primer efficiency you expect: "))
                print efficiency
                break
            except ValueError: 
                print "The value type you entered was invalid. Try again!"
    else:
        sum = 0.0 #ensure they are float type
        count = 0.0 
        for CtDifference in standardCtDifferences:
            calculatedEfficiency = math.pow(10, math.log(float(dilutionFactor), 10.0) / float(standardCtDifferences[CtDifference]))
            sum += calculatedEfficiency
            count += 1
        efficiency = float(sum / count)
    if efficiency > 2.0:
        efficiency = 2.0 #set 2 as maximum primer efficiency
    
    print "Here are the samples you can choose from: " 
    for sample in sampleDictionary.keys():
        startingIndex = sample.index(" ") + 1
        print sample[startingIndex:]
    baselineSample = raw_input("Which sample should act as the baseline? Enter it's name exactly as written above: ")
    baselineSample = primer + " " + baselineSample 
    while baselineSample not in sampleDictionary.keys() or sampleDictionary[baselineSample][0] == "NA": #reprompt if sample is invalid
        baselineSample = primer + " " + raw_input("The name you typed in was either incorrect OR the sample you chose has a 'NA' Ct. Please enter the baseline sample name again: ")
    
    #calculate difference in Ct from baseline AND amount relative to baseline
    baselineCt = sampleDictionary[baselineSample][0]
    relativeToBaseline = {}
    for sample in sampleDictionary.keys():
        if sample != baselineSample: #don't subtract baseline from baseline
            sampleCt = sampleDictionary[sample][0]  
            key = sample
            if sampleCt == "NA":
                relativeToBaseline[key] = ["NA", "NA"]
            else:    
                difference = float(baselineCt) - float(sampleCt)
                amountRelative = math.pow(float(efficiency), difference)
                relativeToBaseline[key] = [round(difference, 6), round(amountRelative, 6)]
        else: #sample == baselineSample
            relativeToBaseline[baselineSample] = ["BASELINE", "BASELINE"]
    
    #Print Out Results in a text file
    F2.write("Primer: " + primer + "\n")
    F2.write("Primer Efficiency: " + str(efficiency))
    F2.write("\n") #skip a line
    F2.write("Sample Name" + "\t" + "Median Ct Value" + "\t" + "Ct Relative to " + baselineSample + "\t" + "Amount Relative to " + baselineSample + "\n")
    for sample in sampleDictionary.keys():
        F2.write(sample + "\t" + str(sampleDictionary[sample][0]) + "\t" + str(relativeToBaseline[sample][0]) + "\t" + str(relativeToBaseline[sample][1]))
        F2.write("\n") #skip a line between each sample 
    return totalSamples
 
startIndex = 0
F2 = open("output.txt", "w")  
for primer in uniquePrimers:
    totalSamples = primerAnalysis(masterSampleDictionary, startIndex, F2)    
    startIndex += totalSamples  
    F2.write("\n") #skip a line between each primer
F2.close()
print "The qPCR analysis is done! Check 'output.txt' for the results."


        
        
    

    


    
