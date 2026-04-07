from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os

# ---------------- USER SETTINGS ----------------
specimenList = ['C3D10_single_Abq']
directory = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/single_C3D10_Abq/'
instanceName = 'PART-1-1'
stepName = 'Step-1'
frequency = 1
SDV_FIELD = 'SDV8'
# -----------------------------------------------


for specimen in specimenList:

    odbPath = os.path.join(directory, specimen + '.odb')
    vtkBase = os.path.join(directory, 'vtk_' + specimen, specimen)

    if not os.path.exists(os.path.dirname(vtkBase)):
        os.makedirs(os.path.dirname(vtkBase))

    print('Opening ODB:', odbPath)
    myOdb = session.openOdb(name=odbPath)

    frames = myOdb.steps[stepName].frames
    numFrames = len(frames)

    myInstance = myOdb.rootAssembly.instances[instanceName]
    nodes = myInstance.nodes
    elements = myInstance.elements

    numNodes = len(nodes)
    numElements = len(elements)

    # ---------------------------------------------------------
    # Build node map and initial coordinates from ODB
    # ---------------------------------------------------------
    InitialCoords = []
    nodeid = {}  # Abaqus label -> VTK index

    for vtk_index, node in enumerate(nodes):
        nodeid[node.label] = vtk_index
        InitialCoords.append(node.coordinates)

    # ---------------------------------------------------------
    # Build element connectivity from ODB
    # ---------------------------------------------------------
    elementConnectivity = []
    elementLabels = []

    for el in elements:
        conn = [nodeid[n] for n in el.connectivity]
        elementConnectivity.append(conn)
        elementLabels.append(el.label)

    # ---------------------------------------------------------
    # Loop over frames
    # ---------------------------------------------------------
    for fr in range(0, numFrames, frequency):

        print('Processing frame:', fr)

        vtkFile = open(vtkBase + str(fr) + '.vtk', 'w')

        vtkFile.write('# vtk DataFile Version 2.0\n')
        vtkFile.write('Abaqus ODB Export\n')
        vtkFile.write('ASCII\n')
        vtkFile.write('DATASET UNSTRUCTURED_GRID\n')

        # -----------------------------------------------------
        # Displacement field mapped by node label
        # -----------------------------------------------------
        if fr > 0:
            uField = frames[fr].fieldOutputs['U'].getSubset(region=myInstance).values
            uDict = {v.nodeLabel: v.data for v in uField}

        # -----------------------------------------------------
        # Write nodes
        # -----------------------------------------------------
        vtkFile.write('POINTS %i float\n' % numNodes)

        for node in nodes:
            x0, y0, z0 = node.coordinates

            if fr == 0:
                x, y, z = x0, y0, z0
            else:
                ux, uy, uz = uDict[node.label]
                x = x0 + ux
                y = y0 + uy
                z = z0 + uz

            vtkFile.write('%f %f %f\n' % (x, y, z))

        # -----------------------------------------------------
        # Write elements (C3D10 -> VTK_QUADRATIC_TETRA = 24)
        # -----------------------------------------------------
        vtkFile.write('CELLS %i %i\n' % (numElements, 11 * numElements))

        for conn in elementConnectivity:
            vtkFile.write('10 ' + ' '.join(str(n) for n in conn) + '\n')

        vtkFile.write('CELL_TYPES %i\n' % numElements)
        for _ in range(numElements):
            vtkFile.write('24\n')

        # -----------------------------------------------------
        # SDV8 averaged per element (mapped by element label)
        # -----------------------------------------------------
        vtkFile.write('CELL_DATA %i\n' % numElements)
        vtkFile.write('SCALARS SDV8 float 1\n')
        vtkFile.write('LOOKUP_TABLE default\n')

        if fr == 0:
            for _ in range(numElements):
                vtkFile.write('0.0\n')
        else:
            sdvField = frames[fr].fieldOutputs[SDV_FIELD].getSubset(region=myInstance).values

            elementVals = {}
            for val in sdvField:
                eid = val.elementLabel
                elementVals.setdefault(eid, []).append(val.data)

            for eid in elementLabels:
                if eid in elementVals:
                    avg = sum(elementVals[eid]) / len(elementVals[eid])
                else:
                    avg = 0.0
                vtkFile.write('%f\n' % avg)

        vtkFile.close()

    myOdb.close()
