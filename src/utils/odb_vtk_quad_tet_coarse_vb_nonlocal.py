# Extract ABAQUS ODB into VTK unstructured grid data format
# Based on the code of Hadi Hosseini, modified to read SDV8 from ODB

# Get ABAQUS interface
from abaqus import *
from abaqusConstants import *
from sys import exit, stdout, argv
from odbAccess import *
from math import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import os

specimenList = ['C3D10_single_Abq']
frequency = 1

DISPLACEMENT = dict()
InitialCoords = dict()

for specimen in specimenList: 
    directory = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/single_C3D10_Abq/'
    odbPath = directory + specimen + '.odb'
    inpPath = directory + specimen + '.inp'
    vtkFile1 = directory + '/vtk_' + specimen + '/' + specimen  

    ##########################################################################################
    instanceName = 'PART-1-1'
    stepName = 'Step-1'
   
    # Open the odb
    myOdb = session.openOdb(name=odbPath)

    # Get the frame repository for the step, find number of frames
    frames = myOdb.steps[stepName].frames
    numFrames = len(frames)
   
    numFreq = numFrames / int(frequency)

    # Isolate the instance, get the number of nodes and elements
    myInstance = myOdb.rootAssembly.instances[instanceName]
    numNodes = len(myInstance.nodes)
    numElements = len(myInstance.elements)

    # Isolate the displacement field
    for fr in range(frequency, numFrames, frequency):
        DISPLACEMENT[fr] = frames[fr].fieldOutputs['U'].getSubset(region=myInstance).values

    ##########################################################################################
    # Get the element connectivity
    elementConnectivity = []
    InitialCoords = []
    nodeid = {}

    with open(inpPath, 'r') as f:
        lines = f.readlines()

    reading_nodes = False
    vtk_index = 0
    for line in lines:
        line = line.strip()
        if line.startswith('*NODE,'):
            reading_nodes = True
            continue
        if reading_nodes and (line.startswith('*') or line.startswith('**')):
            reading_nodes = False
            continue
        if reading_nodes:
            parts = line.split(',')
            abaqus_id = int(parts[0])
            coords = [float(parts[1]), float(parts[2]), float(parts[3])]
            nodeid[abaqus_id] = vtk_index
            InitialCoords.append(coords)  # order matches VTK indices
            vtk_index += 1
      
    elementConnectivity = []

    reading_elements = False
    for line in lines:
        line = line.strip()
        if line.startswith('*ELEMENT,'):
            reading_elements = True
            continue
        if reading_elements and (line.startswith('*') or line.startswith('**')):
            reading_elements = False
            continue
        if reading_elements:
            parts = line.split(',')
            # skip element ID (parts[0]), map all remaining node IDs
            connectivity = [nodeid[int(nid)] for nid in parts[1:]]
            elementConnectivity.append(connectivity)
   
    ##########################################################################################
    # Extract SDV8 values directly from the ODB
    SDV_FIELD = 'SDV8'  # replace with actual SDV name in ODB
    sdvDict = {}  # key: frame number, value: list of element averaged SDV8

    for fr in range(frequency, numFrames, frequency):
        sdvDict[fr] = []

        # Access the SDV field for the current frame
        sdvField = frames[fr].fieldOutputs[SDV_FIELD].getSubset(region=myInstance).values

        # Accumulate SDV8 per element
        elementVals = {}  # key: elementLabel, value: list of SDV8 per gauss point
        for val in sdvField:
            eid = val.elementLabel
            if eid not in elementVals:
                elementVals[eid] = []
            elementVals[eid].append(val.data)  # val.data is SDV8 at Gauss point

        # Average SDV8 per element
        for el in range(1, numElements + 1):
            if el in elementVals:
                avg_sdv = sum(elementVals[el]) / len(elementVals[el])
            else:
                avg_sdv = 0.0
            sdvDict[fr].append(avg_sdv)

    ##########################################################################################
    for fr in range(0, numFrames, frequency):
        # Open the output vtk file and write the header
        vtk_dir = os.path.dirname(vtkFile1)
        if not os.path.exists(vtk_dir):
            os.makedirs(vtk_dir)
        vtkFile = open(vtkFile1 + str(fr) + '.vtk', 'w')

        vtkFile.write('# vtk DataFile Version 2.0\n')
        vtkFile.write('Reconstructed Lagrangian Field Data\n')
        vtkFile.write('ASCII\n')
        vtkFile.write('DATASET UNSTRUCTURED_GRID\n')
        
        # Print out all the coordinates to the vtk file
        vtkFile.write('POINTS %i float\n' % (numNodes))
         
        if fr == 0:
            # write the initial coordinates
            for nd in range(0, numNodes):
                x = InitialCoords[nd][0]
                y = InitialCoords[nd][1]
                z = InitialCoords[nd][2]
                vtkFile.write('%f\t%f\t%f\n' % (x, y, z))        
        else:          
            # Add displacements to the initial coordinates
            for nd in range(0, numNodes):
                x = InitialCoords[nd][0] + DISPLACEMENT[fr][nd].data[0]
                y = InitialCoords[nd][1] + DISPLACEMENT[fr][nd].data[1]
                z = InitialCoords[nd][2] + DISPLACEMENT[fr][nd].data[2]
                vtkFile.write('%f\t%f\t%f\n' % (x, y, z))

        # Print out all the elements to the vtk file
        vtkFile.write('CELLS %i %i\n' % (numElements, 11 * numElements))
     
        for el in range(0, numElements):
            n1, n2, n3, n4, n5, n6, n7, n8, n9, n10 = elementConnectivity[el]
            vtkFile.write('10\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n' % 
                          (n1, n2, n3, n4, n5, n6, n7, n8, n9, n10))

        # Print out all the cell types to the vtk file
        vtkFile.write('CELL_TYPES %i\n' % (numElements))
        for el in range(0, numElements):
            vtkFile.write('24\n')

        vtkFile.write('POINT_DATA %i\n' % (numNodes))
        vtkFile.write('CELL_DATA %i\n' % (numElements))
        
        # Print out the scalar averaged SDV8 field to the vtk file
        vtkFile.write('SCALARS SDV8 float 1\n')   
        vtkFile.write('LOOKUP_TABLE default\n')
        
        for el in range(numElements):
            if fr == 0:
                vtkFile.write('%f\n' % (0.0))
            else:
                vtkFile.write('%f\n' % (sdvDict[fr][el]))

        vtkFile.close()        

    myOdb.close()

    '''
    abaqus viewer noGUI=odb_vtk_quad_tet_coarse_vb_nonlocal.py
    '''