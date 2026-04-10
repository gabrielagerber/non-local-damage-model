from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os

# ---------------- USER INPUT ----------------
specimenList = ['C3D10_cube_NL']
frequency = 1

instanceName = 'PART-1-1'
stepName = 'Step-1'

baseDir = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/cube_C3D10_NL/'
# -------------------------------------------

for specimen in specimenList:

    odbPath = baseDir + specimen + '.odb'
    inpPath = baseDir + specimen + '.inp'
    vtkDir  = baseDir + 'vtk_' + specimen

    if not os.path.exists(vtkDir):
        os.makedirs(vtkDir)

    vtkBase = vtkDir + '/' + specimen

    # ---------------- OPEN ODB ----------------
    odb = session.openOdb(name=odbPath)

    step = odb.steps[stepName]
    frames = step.frames
    numFrames = len(frames)

    instance = odb.rootAssembly.instances[instanceName]

    # ---------------- NODE MAP ----------------
    nodeIdMap = {}
    initialCoords = []

    with open(inpPath, 'r') as f:
        lines = f.readlines()

    reading_nodes = False
    vtk_idx = 0

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
            nid = int(parts[0])
            coords = list(map(float, parts[1:4]))

            nodeIdMap[nid] = vtk_idx
            initialCoords.append(coords)
            vtk_idx += 1

    numNodes = len(initialCoords)

    # ---------------- ELEMENT MAP ----------------
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
            elemLabel = int(parts[0])
            conn = [nodeIdMap[int(n)] for n in parts[1:]]

            elementConnectivity.append((elemLabel, conn))

    numElements = len(elementConnectivity)

    # ---------------- DISPLACEMENTS ----------------
    dispDict = {}

    for fr in range(0, numFrames, frequency):

        frame = frames[fr]
        dispField = frame.fieldOutputs['U'].getSubset(region=instance)

        dmap = {}
        for v in dispField.values:
            dmap[v.nodeLabel] = v.data

        dispDict[fr] = dmap

    # ---------------- SDV8 (AVERAGED OVER IPS) ----------------
    sdv_index = 7  # SDV8
    damageDict = {}

    for fr in range(0, numFrames, frequency):

        frame = frames[fr]
        damageMap = {}

        if 'SDV' in frame.fieldOutputs.keys():

            sdvField = frame.fieldOutputs['SDV'].getSubset(region=instance)

            elemData = {}  # elemLabel -> list of SDV8 from all IPs

            for v in sdvField.values:

                elemLabel = v.elementLabel
                data = v.data

                if len(data) > sdv_index:

                    sdv8 = float(data[sdv_index])

                    if elemLabel not in elemData:
                        elemData[elemLabel] = []

                    elemData[elemLabel].append(sdv8)

            # average over all IPs
            for elemLabel in elemData:
                vals = elemData[elemLabel]
                damageMap[elemLabel] = sum(vals) / len(vals)

        else:
            for el, _ in elementConnectivity:
                damageMap[el] = 0.0

        damageDict[fr] = damageMap

    # ---------------- WRITE VTK ----------------
    for fr in range(0, numFrames, frequency):

        vtkFile = open(vtkBase + str(fr) + '.vtk', 'w')

        vtkFile.write('# vtk DataFile Version 2.0\n')
        vtkFile.write('UEL SDV8 Export (IP averaged)\n')
        vtkFile.write('ASCII\n')
        vtkFile.write('DATASET UNSTRUCTURED_GRID\n')

        # ---- POINTS ----
        vtkFile.write('POINTS %d float\n' % numNodes)

        if fr == 0:
            for i in range(numNodes):
                x, y, z = initialCoords[i]
                vtkFile.write('%f %f %f\n' % (x, y, z))
        else:
            for i in range(numNodes):
                dx, dy, dz = dispDict[fr].get(i+1, (0.0, 0.0, 0.0))
                x, y, z = initialCoords[i]
                vtkFile.write('%f %f %f\n' % (x+dx, y+dy, z+dz))

        # ---- CELLS ----
        vtkFile.write('CELLS %d %d\n' % (numElements, numElements * 11))

        for elemLabel, conn in elementConnectivity:
            vtkFile.write('10 ' + ' '.join(map(str, conn)) + '\n')

        # ---- CELL TYPES (C3D10) ----
        vtkFile.write('CELL_TYPES %d\n' % numElements)

        for _ in range(numElements):
            vtkFile.write('24\n')

        # ---- CELL DATA ----
        vtkFile.write('CELL_DATA %d\n' % numElements)
        vtkFile.write('SCALARS SDV8 float 1\n')
        vtkFile.write('LOOKUP_TABLE default\n')

        for elemLabel, _ in elementConnectivity:
            vtkFile.write('%f\n' % damageDict[fr].get(elemLabel, 0.0))

        vtkFile.close()

    odb.close()

print("Export completed successfully.")

'''
/opt/Abaqus/611/Commands/abq6113 viewer noGUI=odb_vtk_quad_tet_coarse_vb_nonlocal_orig.py
'''