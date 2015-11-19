#Chris Dale
#Bioinformatics Final Project

#Importing for ability to print date and time in file output names
from time import localtime, strftime
#Importing for ability to check if file exists
import os.path

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~     General Methods     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

def printGoal():
    """Prints the purpose of project"""
    print("Chris Dale - Bioinformatics Final Project\n")
    print("The identification of ORF’s from genomic sequences was achieved alone,")
    print("     and with a group Motif identification was added to this program's")
    print("     repertoire of function. \n")
    print("This program can: identify motifs of varying length in whole genomes,")
    print("     identify motifs of varying length in ORFs, and identify motifs of varying")
    print("     length in the 100 bases before an ORF using:")
    print("MoEx - Motif Examiner - returns the actual counts of the occurance ")
    print("     of differnent motifs within the given genome")
    print("ORFMF - Open Reading Frame Motif Finder - returns a count of the ")
    print("     occuances of those motifs within all of the open reading frames ")
    print("     that were found to have a length of over 500 base pairs.")
    print("BSMF - Binding Site Motif Finder - returns a count of the ")
    print("     occuances of the motifs located within the 100 bases before ")
    print("     the start of an open reading frame that was found to have ")
    print("     a length of over 500 base pairs. \n")
    print("I hypothesize that the resulting sequences will probobly represent")
    print("     binding sites of different regulatory protiens. A future program \n")
    print("     could automatically access the online databases and cross reference")
    print("     the sequences with known binding site motifs. \n")
    
def mainWelcome():
    """prints the options one has while using program"""
    print("Options:")
    print("\"MoEx\" for Motif Examiner (full genome), ")
    print("\"ORFMF\" for Open Reading Frame Motif Finder, ")
    print("\"BSMF\" for Binding Site Motif Finder or \"quit\" ")

def openGeneFile(nameFile):
    """generic function for producing string and then list,
    returns a list of all of the characters in the input file"""
    theSequence=open(nameFile,"r",encoding="utf-8")
    theSeqString=theSequence.read()  # read the entire file
    theSequence.close()
    theUpSeqString=theSeqString.upper()
    theSeqList=list(theUpSeqString)
    for i in range(len(theSeqList)-1,-1,-1):
        if theSeqList[i].isalpha()==False or theSeqList[i]=="\t" or theSeqList[i]=="\n":
            del theSeqList[i]
    return(theSeqList)

def giveReverseCompliment(inSequence):
    """returns the reverse compliment of the entire sequence"""
    reverseComp=[]
    for base in range(len(inSequence)):
        if inSequence[base]=="A":
            reverseComp+="T"
        elif inSequence[base]=="T":
            reverseComp+="A"
        elif inSequence[base]=="C":
            reverseComp+="G"
        elif inSequence[base]=="G":
            reverseComp+="C"
    return(reverseComp)


def giveReadFrames(inSequence):
    """returns 3 forward reading frames"""
    leng=len(inSequence)
    seq1=inSequence[0:leng-(leng%3)]
    seq2=inSequence[1:((leng-1)-((leng-1)%3)+1)]
    seq3=inSequence[2:((leng-2)-((leng-2)%3)+2)]
    return(seq1,seq2,seq3)

def isStop(inSequence):
    """end sequences are TGA, TAA and TAG; returns true when given as input"""
    if inSequence==["T","A","A"] or inSequence==["T","A","G"] or inSequence==["T","G","A"]:
        return True
    else:
        return False

def containsMiddleStop(inSequence):
    """Checks for the presence of a stop codon in the middle of the sequence"""
    for tripORF in range(len(inSequence)):
        if(tripORF%3==0 and isStop(inSequence[tripORF:tripORF+3])==True):
            return True
    else:
        return False

def dictIdent(allSeqList,leng,scope):
    """ takes in a list of lists, each sublist is a reading frame, the scope defines
    how the reading frame is broken up into individual ORFs (open reading frames),
    the leng defines the length of ORF being looked for """
    #sets up library and outStr -> string that gets printed to a file in the end
    theSeq={}
    outStr=""
    chimp=0
    #iterated through the 6 sequences
    print("There was not a valkyr in sight...")
    for genome in allSeqList:
        listofbp=[] #empty list
        maximum=0
        #Determines how the information is processed
        if(scope=="BSMF"):
            listORFs=givePreORFs(genome)
        elif(scope=="ORFMF"):
            listORFs=giveORFs(genome)
        elif(scope=="MoEx"):
            listORFs=[genome]
        chimp+=1
        print("And then I found %d %s ..." % (chimp, ['valkyr','valkyrie'][chimp!=1])) #should find 6 nimbi
        for eachORF in listORFs: #gets the ORF's for each reading frame
            for i in range(0,len(genome)-(leng-1),3):
                #Iterating accross ORF in reading frame aviods repeats
                if i%30000==0 and i!=0: #Size of the chunking of the genome
                    keys = [x for x,y in theSeq.items() if y > 1] #add back in the dict if key.value>1
                    for key in keys: #creates a list of keys and values
                        listofbp.append(key)
                        listofbp.append(theSeq[key])
                    theSeq={} #reset dictionary
                    for t in range(0,len(listofbp),2):
                        theSeq[listofbp[t]]=listofbp[t+1]
                bp=genome[i:i+leng]             #sets the motif being looked for
                bpJ="".join(bp)                 #turns motif into string
                if bpJ in theSeq:               #looks for motif string in library
                    theSeq[bpJ]=theSeq[bpJ]+1   #if motif is in lib, adds to value
                else:                           #if motif isn't in lib
                    theSeq[bpJ]=1               #then it adds it to lib with a value of one
    maxxx = max(theSeq.values())
    for h in range(maxxx,1,-1):
        keyss = [x for x,y in theSeq.items() if y==h]
        if(len(keyss)!=0):
            outStr+=str(len(keyss))+" Sequences occurred "+str(h)+" times: "+str(keyss)+" \n \n"
    return outStr

def welcomeGen(thisFile):
    """This is the general welcome function for the motif finding methods"""
    fileEnd=".txt"
    outTime=strftime("%m-%d-%Y_%H-h-%M", localtime())
    fileInA="VCgenomeClean"
    fileInB="HPV11alpha_GUMC-AJ"
    fileInC="lyme_disease"
    print("To choose a sequence, enter either:")
    print(" \"inA\" for Vibrio cholera")
    print(" or \"inB\" for Human Papyloma Virus 11-alpha")
    print(" or \"inC\" for lyme disease")
    print(" or a file name (without .txt) or \"quit\"")
    fil=input("What sequence would you like to process? \n")
    fileFinal=fil+"-Error-Quit_at_time-"+outTime+fileEnd
    leng=0
    #file should never have this name, error if found b/c it means it made a file after a quit
    if(fil=="quit"): #if quitting, returns here
        return fil, leng, fileFinal
    numIn=input("What is the motif length? Type \"d\" for default (15) or enter number. \n")
    if(numIn=="default" or numIn=="d"):
        leng=15
    else:
        while(numIn.isdigit()==False):
            numIn=input("Error factor must be a number, try again. \n")
        leng=int(numIn)
    print("To choose an output, enter either:")
    print(" \"d\" for default ")
    print(" or a file name (without .txt)")
    filOut=input("What would you like to do? \n")
    if(filOut=="d"):
        if(fil=="inA"):
            fileFinal=thisFile+"_Output_Sequence-length="+str(leng)+"-fromFile_"+fileInA+"-"+outTime+fileEnd
        elif(fil=="inB"):
            fileFinal=thisFile+"_Output_Sequence-length="+str(leng)+"-fromFile_"+fileInB+"-"+outTime+fileEnd
        elif(fil=="inC"):
            fileFinal=thisFile+"_Output_Sequence-length="+str(leng)+"-fromFile_"+fileInC+"-"+outTime+fileEnd
        else:
            fileFinal=thisFile+"_Output_Sequence-length="+str(leng)+"-fromFile_"+fil+"-"+outTime+fileEnd
    else:
        fileFinal=thisFile+"_Output_Sequence-length="+str(leng)+"-fromFile_"+filOut+"-"+outTime+fileEnd

    if(fil=="inA"):
        return fileInA+fileEnd, leng, fileFinal
    elif(fil=="inB"):
        return fileInB+fileEnd, leng, fileFinal
    elif(fil=="inC"):
        return fileInC+fileEnd, leng, fileFinal
    else:
        return fil+fileEnd, leng, fileFinal

def runGen(theMethod):
    keepRunning=True
    quitting=False
    while(keepRunning):
        fileNameClean, leng, fileNameFinal=welcomeGen(theMethod)
        while(not os.path.isfile(fileNameClean)):
            if(fileNameClean=="quit"):
                quitting=True
                break
            else:
                print("Sorry, file not found!")
                fileNameClean, leng, fileNameFinal=welcomeGen(theMethod)
        if(quitting):
            break
        seqF=openGeneFile(fileNameClean)
        seqR=giveReverseCompliment(seqF)
        seqF1,seqF2,seqF3=giveReadFrames(seqF)
        seqR1,seqR2,seqR3=giveReadFrames(seqR)
        fileClean=open(fileNameFinal, "w")
        allSeqList=[seqF1,seqF2,seqF3,seqR1,seqR2,seqR3]
        """ takes in the each of the 6 sequence's ORFs and finds the
        DNA sequences that are repeated most often immediatly before
        the ORFs
        """
        outStr=dictIdent(allSeqList,leng,theMethod)
        print("Data written to "+fileNameFinal)
        fileClean.write(outStr)
        check=input("Should I Keep Running? type \"n\" to exit"+theMethod)
        if(check=="n" or check=="no" or check=="N" or check=="NO" or check=="No"):
            keepRunning=False

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      MoEx Methods      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""  

def doMoEx():
    print("Welcome to MoEx - The Motif Examiner ")
    print("This program identifies common motifs among the 6 different reading frames. \n")
    runGen("MoEx")

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   ORFMF Methods     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~""" 

def giveORFs(inSequence):
    """returns the number of ORFs in inSequence"""
    outSeqs=[]
    for trip in range(len(inSequence)):
        if trip%3==0 and inSequence[trip:trip+3]==["A","T","G"]:
            counter=0
            iterator=trip
            done=False
            while iterator<len(inSequence) and done!=True:
                if(iterator%3==0 and counter>500 and isStop(inSequence[iterator:iterator+3])):
                    outSeqs.append(inSequence[trip:iterator+3])
                elif(iterator%3==0 and counter<500 and isStop(inSequence[iterator:iterator+3])):
                    done=True
                iterator+=1
                counter+=1
    return outSeqs
    
def doORFMF():
    print("Welcome to BSMF – Binding Site Motif Finder ")
    print("This program identifies Open Reading Frames and looks to the preceeding ")
    print("genetic code to identify commonly repeated motifs, I hypothesize that the ")
    print("resulting sequences will probobly represent binding sites of different ")
    print("regulatory protiens. \n")
    runGen("ORFMF")

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   BSMF Methods     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~""" 

def givePreORFs(inSequence):
    """returns the bases preceeding a start codon for a full length ORF"""
    outSeqs=[]
    for trip in range(len(inSequence)):
        if trip%3==0 and inSequence[trip:trip+3]==["A","T","G"]:
            counter=0
            iterator=trip
            done=False
            while iterator<len(inSequence) and done!=True:
                if(iterator%3==0 and counter>500 and isStop(inSequence[iterator:iterator+3])):
                    if(trip-100<0):
                        outSeqs.append(inSequence[0:trip])
                    else:
                        outSeqs.append(inSequence[trip-100:trip])
                elif(iterator%3==0 and counter<500 and isStop(inSequence[iterator:iterator+3])):
                    done=True
                iterator+=1
                counter+=1
    return outSeqs 

def doBSMF():
    print("Welcome to BSMF – Binding Site Motif Finder ")
    print("This program identifies Open Reading Frames and looks to the preceeding ")
    print("genetic code to identify commonly repeated motifs. \n")
    runGen("BSMF")

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Main  Method     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~""" 

def main():
    print("Welcome to Chris Dale's Final Project")
    printGoal()
    mainWelcome()
    while(True):
        fil=input("What would you like to do?\n")
        if(fil=="AACC"):
            doAACC()
        elif(fil=="ORFMF"):
            doORFMF()
        elif(fil=="BSMF"):
            doBSMF()
        elif(fil=="MoEx"):
            doMoEx()
        elif(fil=="quit"):
            print("Quitting...")
            break
        elif(fil=="goal"):
            printGoal()
        else:
            mainWelcome()
    print("Program Closing, Goodbye!")
    exit
        

main()
