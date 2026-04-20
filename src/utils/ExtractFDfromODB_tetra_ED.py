"""
#########################################################################
File   : abqFuReader.py (former rcuReader.py)
Author : D.H.Pahr-Modified E. Dall'Ara
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To calculate the displacement of the center of the loading plate and the sum of the reaction forces in 

                -odb     odbName
                -rpt     reportFileName
                [-exe]   abaqusCommand    (optional)
                [-step]  stepID           (optional)
                [-frame] frameID          (optional)
                [-inst]  instanceName     (optional)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Arguments 
 -odb   :  Name of the output database.
 
 -rpt   :  Name of the output file.
 
 -exe   :  Abaqus command which should be used to call this script 
           (needed for medtool)

 -step  :  Step which should be read (default read ALL steps). The steps 
           are given within double dots e.g. 2:6:7:10
             
 -frame :  Frame which should be read (default frame 1). The frame are
           given within double dots e.g. 1:3:4:7. A further possibility
           is the option 'all'.

 -inst  :  Instance name 
           Default is 'PART-1-1'

 -help  :  Print usage
#########################################################################
"""

abqImportError = 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try:
  from odbAccess import *
except ImportError:
  abqImportError = 1 

from sys import argv,exit,stdout
#from Numeric import *    # for array
from numpy import *    # for array
from   time    import *

import os


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def rightTrim(input,suffix):
    if (input.find(suffix) == -1):
        input = input + suffix
    return input
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def getData(odbName, rptName, stepID, frameID, instName='PART-1-1'):
    """ Get  CF, RF, and U: given odbName and rptName"""
    # Master node sets 
    #nsets = 'SEB','SWB','NWB','SET','SWT','NWT'
    nsets = 'NODLOADPLATE','ALL_NODE_T'
    #print nsets
    region = "over the entire model"
    # Open the output database
    odb = openOdb(odbName)
    #try :
    #    odb = openOdb(odbName)
    #except  OdbError :
    #    stdout.write("\n **ERROR** -odb: intput file '%s' not found!\n\n" % odbName ); stdout.flush()
    #    stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
    #    exit(1)
    # Open report file 
    rptOS  = open(rptName, 'w')
    rpt3SName = name.replace(".odb","_Lvdts.txt")
    rpt3S = open(rpt3SName, 'w')
    Data = name.replace(".odb","_FD_NEWBC_tetra_Avg.txt") 
    DataOS = open(Data, 'w')
    rptOS.write('############################################\n')
    rptOS.write('### ECHO of  Coord, CF, RF, and U \n')
    rptOS.write('### odbfile : '+ odbName + '\n')
    rptOS.write('############################################\n')
    #print odb.rootAssembly.instances.keys()
    #print odb.rootAssembly.instances.has_key('PART-1-1')
    if not odb.rootAssembly.instances.has_key(instName) :
      stdout.write( "\n **ERROR** Instance %s not found in odb file!" % (instName) )
      stdout.write( "\n           Available: %s \n\n" % odb.rootAssembly.instances.keys() )
      stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
      exit(0)

    myInstance = odb.rootAssembly.instances[instName]
    
    #print myInstance.nodeSets.keys()
    
    # Loop over all master node sets, check nsets 
    for nset in nsets :
        try:
            nsetObj = myInstance.nodeSets[nset]
            region = " Found the nodeset : " + nset;
            print region 
        except KeyError:
            print '** ERROR ** An Instance level NSET named %s does not exist in the '\
                  'output database %s' % (nset, odbName)
            odb.close()
            exit(0)
            
    # Initialize  values 
    maxStep = "_None_"
    maxFrame = -1
    #ConForce = 'CF'
    ReactForce = 'RF'
    Displacement = 'U' 
    noStep = 0
    uList = {}
    xyzList = {}
    cfList = {}
    rfList = {}

    #extract step and frame values to list
    stepIDs=[]
    if stepID :
      stepID2 = stepID.replace(':', ' ')
      sStepIDs =  stepID2.split()
      for sStepID in  sStepIDs :
         stepIDs.append( int(sStepID) )

    frameIDs=[]
    if frameID :
      if frameID != 'all' : 
        frameID2 = frameID.replace(':', ' ')
        sFrameIDs = frameID2.split()
        for sFrameID in  sFrameIDs :
           frameIDs.append( int(sFrameID) )
    
    # Loop over all steps
    counter = 0
    for step in odb.steps.values():
      noStep = noStep + 1
      # check if step ids are given, if yes process only given steps
      # if no process all steps 
      if (stepIDs.count(noStep) > 0) or (stepID==None) :
        # Loop over all frames of step, check if frame ids are given
        # if no process only frame 1
        if frameID == 'all' :
           for frame in step.frames:
              frameIDs.append( frame.frameId )
        for frame in step.frames:           
          if (  (frameIDs.count(frame.frameId) > 0) or ( (frameID==None)  and (frame.frameId == 1) )  ):
            print ' ... Processing Step:', step.name, '   Frame:', frame.frameId
            allFields = frame.fieldOutputs
            
            if ( allFields.has_key(ReactForce) ) :
                allRFSet = allFields[ReactForce]
            else : 
              stdout.write( "\n **ERROR** Variable 'RF' not found in odb file!" )
              stdout.write( "\n           Use '*NODE OUTPUT' keyword in inp file!\n\n" )
              stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
              exit(0)
            if ( allFields.has_key(Displacement) ) :
                allUSet  = allFields[Displacement]
            else : 
              stdout.write( "\n **ERROR** Variable 'U' not found in odb file!" )
              stdout.write( "\n           Use '*NODE OUTPUT' keyword in inp file!\n\n" )
              stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
              exit(0)



            # COORDINATES
            counter = counter + 1 
            for nset in nsets :
               nsetObj = myInstance.nodeSets[nset]
               USet   = allUSet.getSubset(region=nsetObj)
               for UValue in USet.values:
                    curNode  = UValue.nodeLabel
                    if counter == 1 : 
                       # save coordinates only once at the first step 
                       vXYZ=array((myInstance.getNodeFromLabel(curNode).coordinates[0],
                                  myInstance.getNodeFromLabel(curNode).coordinates[1],
                                  myInstance.getNodeFromLabel(curNode).coordinates[2]))   
                       xyzList[nset] = vXYZ
                       #strList = str(xyzList[nset])
                       #wstr    = str('XYZ; ' + nset + '; ' +  strList + '\n')
                       #rptOS.write( wstr )
                       #print 'XYZ: ',nset, xyzList[nset]   

            # RESULTS
            nsetStep = '---:' + step.name + ':Frame-' + str(frame.frameId) 
            wstr    = str('TIME; ' + nsetStep + '; ' +  str(frame.frameValue) + '\n')
            rptOS.write( wstr ) 
	    rpt3S.write( wstr ) 

            for nset in nsets :
	      if nset == 'NODLOADPLATE' :
               nsetObj = myInstance.nodeSets[nset]
               RFSet  = allRFSet.getSubset(region=nsetObj)
               USet   = allUSet.getSubset(region=nsetObj)
               nsetStep = nset + ':' + step.name
               vUAvg  = 0
	       vRFSum = 0
               for UValue in USet.values:
                    curNode  = UValue.nodeLabel

                    # save displacement 
                    vU=array((UValue.data[0],UValue.data[1],UValue.data[2]))
		    uList[nsetStep] = vU
		    vUAvg= vUAvg+(UValue.data[2])/3
		    strList = str(uList[nsetStep])
		    strAvg  = str(vUAvg)
		    #wwstr    = str('ULvdts ; ' + nsetStep + '; ' +  strList + '\n') 
		    wwstr    = str(strList + '\n') 
		    #Avgstr  = str('UAvg   ; ' + nsetStep + '; ' +  strAvg  + '\n')
                    rpt3S.write("%1.5f\t%1.5f\t%1.5f\n" % ( UValue.data[0], UValue.data[1], UValue.data[2]) )
		    #rpt3S.write( wwstr ) 
		    #rptOS.write( Avgstr )                 
                    # curNode  = UValue.nodeLabel
                    # curStep  = step.name
                    # curFrame = frame.frameId

               #for RFValue in RFSet.values:
                    #vRF=array(( RFValue.data[0], RFValue.data[1], RFValue.data[2]))
                    #rfList[nsetStep] = vRF
                    #strList = str(rfList[nsetStep])
                    #wstr    = str('RF  ; ' + nsetStep + '; ' +  strList + '\n') 
                    #rptOS.write( wstr )                          
              else:
               nsetObj = myInstance.nodeSets[nset]
               RFSet  = allRFSet.getSubset(region=nsetObj)
               USet   = allUSet.getSubset(region=nsetObj)
               nsetStep = nset + ':' + step.name

               for RFValue in RFSet.values:
                    vRF=array(( RFValue.data[0], RFValue.data[1], RFValue.data[2]))
		    vRFSum= vRFSum+RFValue.data[2]
                    #rfList[nsetStep] = vRF
                    #strList = str(rfList[nsetStep])
                    #wstr    = str('RF  ; ' + nsetStep + '; ' +  strList + '\n') 
		    strRFSum = str(vRFSum)
                    #wstr      = str('SumRF  ; ' + nsetStep + '; ' +  vRFSumstr + '\n')
                    #rptOS.write( wstr ) 
	    #rptOS.write( Avgstr )
	    #rptOS.write( vRFSumstr )                          
            Datastr  = str('UAvg-RFsum; ' + nsetStep + '; ' +  strAvg + '; ' +  strRFSum  + '\n')
	    rptOS.write(Datastr)
	    Datastr2 = str( strAvg + ' ' +  strRFSum  + '\n')
	    DataOS.write(Datastr2)
	
	
	###convert the file in a table for Mathe
      
      #PointAx = []
      #PointAy = []
      #PointAz = []
      #PointBx = []
      #PointBy = []
      #PointBz = []
      #PointCx = []
      #PointCy = []
      #PointCz = []
      
      #count=0
      

      #for line in open(rpt3SName,'r'):
        
	#if count==3:
	  #count=0       
        #else: pass
	
	###if line[0] == '[' :
	#if line[0] != 'T' :
          #count=count+1
	  ### read this line 
	  
	  #print line
	  
	  ##line = line.replace("]", " ")
	  ##line = line.replace("[   ", "")
	  ##line = line.replace("[  ", "")
	  ##line = line.replace("[ ", "")
	  ##line = line.replace("[", "")
	  ##line = line.replace("                 ", "\t")
	  ##line = line.replace("                ", "\t")
	  ##line = line.replace("               ", "\t")
	  ##line = line.replace("              ", "\t")
	  ##line = line.replace("             ", "\t")
	  ##line = line.replace("            ", "\t")
	  ##line = line.replace("           ", "\t")
	  ##line = line.replace("          ", "\t")
	  ##line = line.replace("         ", "\t")
	  ##line = line.replace("        ", "\t")
	  ##line = line.replace("       ", "\t")
	  ##line = line.replace("      ", "\t")
	  ##line = line.replace("     ", "\t")
	  ##line = line.replace("    ", "\t")
	  ##line = line.replace("   ", "\t")
	  ##line = line.replace("  ", "\t")
          ##line = line.replace(" ", "\t")
	  ##line = line.replace("\n", "\t")
          ##line = line.replace("\r", "")
          ##line = line.replace(",", ".")
	  
	  ##print line
	  
	  ##valueList = line.split("\t")
	  
	  ##if count==1: 
	    ##PointAx.append(float(valueList[0]))
            ##PointAy.append(float(valueList[1]))
	    ##PointAz.append(float(valueList[2]))
	  
	  ##elif count==2: 
	    ##PointBx.append(float(valueList[0]))
            ##PointBy.append(float(valueList[1]))
	    ##PointBz.append(float(valueList[2]))
	  
	  ##elif count==3: 
	    ##PointCx.append(float(valueList[0]))
            ##PointCy.append(float(valueList[1]))
	    ##PointCz.append(float(valueList[2]))
	    ##print count
          ##else :
          ## skip this line 
          ## pass
	
      #PointAx = array(PointAx)
      #PointAy = array(PointAy)
      #PointAz = array(PointAz)
      #PointBx = array(PointBx)
      #PointBy = array(PointBy)
      #PointBz = array(PointBz)
      #PointCx = array(PointCx)
      #PointCy = array(PointCy)
      #PointCz = array(PointCz)
       
       
      #dataIS = open(rpt3SName.replace(".txt","_da1.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointAx[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_da2.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointAy[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_da3.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointAz[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_db1.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointBx[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_db2.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointBy[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_db3.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointBz[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_dc1.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointCx[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_dc2.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointCy[i]))
      
      #dataIS = open(rpt3SName.replace(".txt","_dc3.txt"), 'w')
      #for i in range(len(PointAx)):
        #dataIS.write(" %1.5f \n" % (PointCz[i]))
       
	    
    """ Close the output database before exiting the program """
    odb.close()
    
#==========================================================================
# S T A R T
#    
if __name__ == '__main__':
    
    odbName  = None
    rptName  = None
    stepID   = None
    frameID  = None
    instName = None

    # Type             Param Description               DefaultValue   Optional    Info"
    guiInit = \
   "*modulName        'abqFuReader'                                                        \n"\
   +"*fileEntryIn     -odb   'Input File Name'         testIn.odb     no          odb      \n"\
   +"*fileEntryOut    -rpt   'Output File Name'        testOut.rpt    no          txt:rpt  \n"\
   +"*entry           -exe   'Abaqus Command'          /programs/hks/Commands/abq:python  yes   1 \n"\
   +"*entry           -step  'Analysis Step'           1:2:3:4:5:6    yes         1        \n"\
   +"*entry           -frame 'Analysis Frame'          1:3:4          yes         1        \n"\
   +"*entry           -inst  'Instance Name'           PART-1-1       yes         1        \n"

    argList = argv
    argc = len(argList)
    i=0
    while (i < argc):
        if (argList[i][:4] == "-odb"):
            i += 1
            name = argList[i]
            odbName = rightTrim(name,".odb")
        elif (argList[i][:4] == "-rpt"):
            i += 1
            rptName = argList[i]
        elif (argList[i][:4] == "-ste"):
            i += 1
            stepID  = argList[i]
        elif (argList[i][:4] == "-fra"):
            i += 1
            frameID = argList[i]
        elif (argList[i][:5] == "-inst"):
            i += 1
            instName = argList[i]
        elif (argList[i][:4] == "-gui"):
            stdout.write( '%s' % guiInit); stdout.flush()
            exit(0)
        elif (argList[i][:2] == "-h"):            
            print __doc__
            exit(0)
        i += 1
    
    stdout.write( '\n S T A R T  abqFuReader.py  V_08.04.2008- D. H. Pahr \n\n' )
    stdout.flush()
                    
    if not (odbName):
        stdout.write( __doc__ ); stdout.flush()
        stdout.write( "\n **ERROR** '-odb' file name not given\n\n" ); stdout.flush()
        stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
        exit(1)
    if not (rptName):
        stdout.write( __doc__ ); stdout.flush()
        stdout.write( "\n **ERROR** '-rpt' file name not given\n\n" ); stdout.flush()
        stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
        exit(1)
    if abqImportError: 
        stdout.write( "\n **ERROR** ABAQUS python modules could not be loaded!\n\n" ); stdout.flush()
        stdout.write( '\n E N D E D  with ERRORS \n\n' ); stdout.flush()
        exit(1)      

    #######################################################################

    
    startTime = clock()
    ctime1    = time()
    
    if instName == None : instName = 'PART-1-1'

    getData(odbName,rptName,stepID,frameID,instName)
    endTime = clock()
    ctime2 = time()
    stdout.write( '\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime-startTime, ctime2-ctime1 ) ); stdout.flush()

'''
/opt/Abaqus/611/Commands/abq6113 python ExtractFDfromODB_tetra_ED.py -odb 02121938L1_105_233_QUAEQ_c5_tet_uel.odb -rpt 02121938L1_105_233_QUAEQ_c5_tet_uel_Output.txt -frame all
'''
    ###/programs/hks/Commands/abq681 python ExtractFDfromODB_tetra_ED.py -odb 02121938L1.odb -rpt 02121938L1_Output.txt -frame all
    ###/programs/hks/Commands/abq681 python ExtractFDfromODB_tetra_ED.py -odb $filename.odb -rpt $filename_Output.txt -frame all
