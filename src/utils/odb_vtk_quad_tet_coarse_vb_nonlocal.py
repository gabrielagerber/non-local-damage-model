# Extract ABAQUS ODB into VTK unstructured grid data format

# Get ABAQUS interface
from abaqus import *
from abaqusConstants import *
from sys import exit, stdout,argv
from odbAccess import *
from math import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import os

loadingElement=1
loadingNodes=4

specimenList=['C3D10_single']

frequency = 1

DISPLACEMENT=dict()
InitialCoords=dict()

for specimen in specimenList: 
   directory = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/00_Original_Files/HadiHosseini_March2026/'
   odbPath = directory+specimen+'.odb'
   CtxInpPath = directory+specimen+'.inp'
   TrabInpPath = directory+specimen+'.inp'
   datPath = directory+specimen+'.dat'
   vtkFile1 = directory+'/vtk_'+specimen+'/'+specimen+'_coarse_nonlocal_'   

##########################################################################################
   instanceName = 'PART-1-1'
   stepName = 'Step-1'
   
   # Open the odb
   myOdb = session.openOdb(name=odbPath)
   #odb=openOdb(path = odbPath, readOnly=True)

   # Get the frame repository for the step, find number of frames
   frames = myOdb.steps[stepName].frames
   numFrames = len(frames)
   
   numFreq = numFrames/int(frequency)

   # Isolate the instance, get the number of nodes and elements
   myInstance = myOdb.rootAssembly.instances[instanceName]
   numNodes1 = len(myInstance.nodes)
   numNodes=numNodes1-loadingNodes
   numElements1 = len(myInstance.elements)
   numElements=numElements1-loadingElement
   print ('numNodes = ', numNodes)
   print ('numElements = ', numElements)

   # Isolate the displacement field
   for fr in range(frequency,numFrames,frequency):
     DISPLACEMENT[fr]=dict()  
     DISPLACEMENT[fr]=frames[fr].fieldOutputs['U'].getSubset(region=myInstance).values
        
   # Get the initial nodal coordinates
   for nd in range(0, numNodes):
     coords = myInstance.nodes[nd].coordinates
     InitialCoords[nd] = []
     InitialCoords[nd].append(coords[0])
     InitialCoords[nd].append(coords[1])
     InitialCoords[nd].append(coords[2])
     #print 'nd = ', nd, ' Initial Coordinates = ', InitialCoords[nd][0], InitialCoords[nd][1], InitialCoords[nd][2]
     #print 'nd = ', nd, ' DISPLACEMENT[100]=', DISPLACEMENT[100][nd].data[0],DISPLACEMENT[100][nd].data[1],DISPLACEMENT[100][nd].data[2] 
##########################################################################################
   # Get the element connectivity
   elementConnectivity = []
   nodeid = {}
   CtxInpFile = open(CtxInpPath,'r')
   Ctxlinesinp = CtxInpFile.readlines()
   TrabInpFile = open(TrabInpPath,'r')
   Trablinesinp = TrabInpFile.readlines()   
   i=0
   switch=0
   switch3=0
   while i <len(Ctxlinesinp):

      if Ctxlinesinp[i].find('*NODE')>-1 and Ctxlinesinp[i].find('*NSET')>-1: 
         switch3=1  
         ii=i   
      if switch3 == 1 and Ctxlinesinp[i].find('*NODE') == -1 and Ctxlinesinp[i].find('**') == -1:
         node =  Ctxlinesinp[i].split(',')
         nodeid[int(node[0])]=i-ii-1
      if Ctxlinesinp[i].find('***') >-1: 
         switch3 =0
         print ('dict length = ', len(nodeid))
      
      if Ctxlinesinp[i].find('*ELEMENT,')>-1 and Ctxlinesinp[i].find('TYPE,')>-1: switch=1
      if switch==1: 
         con = []
         element = Ctxlinesinp[i+1].split(',')
           
         elementConnectivity.append((nodeid[int(element[1])],nodeid[int(element[2])],nodeid[int(element[3])],\
         nodeid[int(element[4])],nodeid[int(element[5])],nodeid[int(element[6])],nodeid[int(element[7])],\
         nodeid[int(element[8])],nodeid[int(element[9])],nodeid[int(element[10])]))
      switch=0
      i=i+1      
   CtxInpFile.close()  
   
##########################################################################################   
   # Get the damage values
   datFile = open(datPath,'r')
   linesdat = datFile.readlines()
      
   damageDict = {} # key: time increment id and element id, value: averaged damage over 4 Gauss points

   elementDict = {}

   i=1
   l=i
   for inc in range (frequency,numFrames,frequency):
 
       switch2=0
       switch1=0
       while l <len(linesdat) and switch2==0:
         if linesdat[i].find('INCREMENT     '+str(inc)+' SUMMARY')>-1 or linesdat[i].find('INCREMENT    '+\
         str(inc)+' SUMMARY')>-1 or linesdat[i].find('INCREMENT   '+str(inc)+' SUMMARY')>-1 or linesdat[i].find('INCREMENT  '+str(inc)+' SUMMARY')>-1:
            l=i
            switch2=1
            switch1=1
            j=i+1
            k=0
            while j <len(linesdat) and switch1==1:
               if linesdat[j].find('SUMMARY')>-1: 
                  switch1=0
               if linesdat[j].find('ELEMENT  PT FOOT')>-1 and switch1==1:
                  SDV = linesdat[j+3].split()
                  damageGP = [float(x) for x in SDV]
                  k=k+1
                  if not damageGP: 
                     damage = 0.0
                     elementDict={k:damage}
                     damageDict.setdefault(inc,[]).append(elementDict)
                  else:
                     damage = (damageGP[2]+damageGP[3]+damageGP[4]+damageGP[5])/4
                     elementDict={k:damage}
                     damageDict.setdefault(inc,[]).append(elementDict)
               j=j+1       
         i=i+1

   datFile.close()   
##########################################################################################
   for fr in range(0,numFrames,frequency):
        # Open the output vtk file and write the header
        print ('i am here 1')
        vtk_dir = os.path.dirname(vtkFile1)  # extracts the folder path
        if not os.path.exists(vtk_dir):
            os.makedirs(vtk_dir)
        vtkFile = open(vtkFile1+str(fr)+'.vtk', 'w')
        print ('i am here 2')
        vtkFile.write('# vtk DataFile Version 2.0\n')
        print ('i am here 3')
        vtkFile.write('Reconstructed Lagrangian Field Data\n')
        vtkFile.write('ASCII\n')
        vtkFile.write('DATASET UNSTRUCTURED_GRID\n')
        
        # Print out all the coordinates to the vtk file
        vtkFile.write('POINTS %i float\n' % (numNodes))
         
        if (fr==0):
           # write the initial coordinates
           for nd in range(0, numNodes):
               x = InitialCoords[nd][0]
               y = InitialCoords[nd][1]
               z = InitialCoords[nd][2]
               vtkFile.write('%f\t%f\t%f\n' % (x, y, z))        
        else:          
           # Add displacements to the initial coordinates
           for nd in range(0, numNodes):
               x = InitialCoords[nd][0]+DISPLACEMENT[fr][nd].data[0]
               y = InitialCoords[nd][1]+DISPLACEMENT[fr][nd].data[1]
               z = InitialCoords[nd][2]+DISPLACEMENT[fr][nd].data[2]
               vtkFile.write('%f\t%f\t%f\n' % (x, y, z))

        # Print out all the elements to the vtk file
        vtkFile.write('CELLS %i %i\n' % (numElements, 11*numElements))
     
        for el in range(0,numElements):
            n1  = elementConnectivity[el][0]
            n2  = elementConnectivity[el][1]
            n3  = elementConnectivity[el][2]
            n4  = elementConnectivity[el][3]
            n5  = elementConnectivity[el][4]
            n6  = elementConnectivity[el][5]
            n7  = elementConnectivity[el][6]
            n8  = elementConnectivity[el][7]
            n9  = elementConnectivity[el][8]
            n10 = elementConnectivity[el][9]
            vtkFile.write('10\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n' % (n1, n2, n3, n4, n5, n6, n7, n8, n9, n10))

        # Print out all the cell types to the vtk file
        vtkFile.write('CELL_TYPES %i\n' % (numElements))
        for el in range(0,numElements):
            vtkFile.write('24\n')

        vtkFile.write('POINT_DATA %i\n' % (numNodes))

        vtkFile.write('CELL_DATA %i\n' % (numElements))
        
        # Print out the scalar averaged damage field to the vtk file
        
        vtkFile.write('SCALARS damage float 1\n')   
        vtkFile.write('LOOKUP_TABLE default\n')
        
        for el in range(1,numElements+1):
            if fr==0:
               vtkFile.write('%f\n' % (0.0))
            else:
               vtkFile.write('%f\n' % (damageDict[fr][el-1][el]))

        vtkFile.close()        

   myOdb.close()

'''
/opt/Abaqus/611/Commands/abq6113 viewer noGUI=odb_vtk_quad_tet_coarse_vb_nonlocal.py
'''
