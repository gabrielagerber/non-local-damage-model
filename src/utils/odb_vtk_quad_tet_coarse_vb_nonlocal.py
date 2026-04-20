'''
This script is based on the work of Hadi Hosseini.
It can be used to visualise the Abaqus simulation outcomes using UEL.

CAVE:
All variables expressed at integration points (S,LE,SDV) are averaged across the element.
'''

from abaqus import *
from abaqusConstants import *
from odbAccess import *
import os

specimenList = ['C3D10_single']
frequency = 1  # frame frequency

directory = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/single_C3D10_L/'

DISPLACEMENT = dict()  # nodal displacements
InitialCoords = dict()

# -------------------------------
# Helper function: average element field (scalar or tensor)
def average_element_field(fieldValues, numElements):
    """
    Returns a list of length numElements, each value is the average
    over all Gauss points for that element.
    Supports scalar, vector (tuple/list), and tensor (tuple/list of 6)
    """
    elementVals = {}
    for val in fieldValues:
        eid = val.elementLabel
        if eid not in elementVals:
            elementVals[eid] = []
        if isinstance(val.data, float):
            elementVals[eid].append(val.data)
        else:
            elementVals[eid].append(val.data)

    avg_list = []
    for el in range(1, numElements + 1):
        if el in elementVals:
            vals = elementVals[el]
            n = float(len(vals))
            # compute component-wise average
            if isinstance(vals[0], float):
                avg_val = sum(vals) / n
            else:
                # tuple/list
                avg_val = [sum([v[i] for v in vals]) / n for i in range(len(vals[0]))]
        else:
            # default zero
            if isinstance(fieldValues[0].data, float):
                avg_val = 0.0
            else:
                avg_val = [0.0]*len(fieldValues[0].data)
        avg_list.append(avg_val)
    return avg_list

# -------------------------------
# Detect if field is nodal
def is_nodal_field(fieldValues):
    try:
        return hasattr(fieldValues[0], 'nodeLabel') and fieldValues[0].nodeLabel is not None
    except:
        return False

# -------------------------------
for specimen in specimenList:
    odbPath = directory + specimen + '.odb'
    inpPath = directory + specimen + '.inp'
    vtkFile1 = directory + '/vtk_' + specimen + '/' + specimen

    instanceName = 'PART-1-1'
    stepName = 'Step-1'

    # Open ODB
    myOdb = session.openOdb(name=odbPath)
    frames = myOdb.steps[stepName].frames
    numFrames = len(frames)

    myInstance = myOdb.rootAssembly.instances[instanceName]
    numNodes = len(myInstance.nodes)
    numElements = len(myInstance.elements)

    # -------------------------------
    # Extract nodal displacements
    for fr in range(frequency, numFrames, frequency):
        if 'U' in frames[fr].fieldOutputs.keys():
            DISPLACEMENT[fr] = frames[fr].fieldOutputs['U'].getSubset(region=myInstance).values

    # -------------------------------
    # Read nodes and elements from .inp
    nodeid = {}
    InitialCoords = []
    vtk_index = 0
    with open(inpPath, 'r') as f:
        lines = f.readlines()

    reading_nodes = False
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
            InitialCoords.append(coords)
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
            connectivity = [nodeid[int(nid)] for nid in parts[1:]]
            elementConnectivity.append(connectivity)

    # -------------------------------
    # Extract all field variables
    fieldDict = {}  # frame -> field_name -> list of values
    for fr in range(frequency, numFrames, frequency):
        fieldDict[fr] = {}
        frame_fields = frames[fr].fieldOutputs.keys()

        for fld in frame_fields:
            fieldValues = frames[fr].fieldOutputs[fld].getSubset(region=myInstance).values
            vals_per_element_or_node = []

            if is_nodal_field(fieldValues):
                # Nodal field: store vector or scalar
                for val in fieldValues:
                    vals_per_element_or_node.append(val.data)
                while len(vals_per_element_or_node) < numNodes:
                    vals_per_element_or_node.append(0.0)
            else:
                # Element-based field: average over gauss points
                vals_per_element_or_node = average_element_field(fieldValues, numElements)
                while len(vals_per_element_or_node) < numElements:
                    if isinstance(fieldValues[0].data, float):
                        vals_per_element_or_node.append(0.0)
                    else:
                        vals_per_element_or_node.append([0.0]*len(fieldValues[0].data))

            fieldDict[fr][fld] = vals_per_element_or_node

    # -------------------------------
    # Write VTK files per frame
    for fr in range(0, numFrames, frequency):
        vtk_dir = os.path.dirname(vtkFile1)
        if not os.path.exists(vtk_dir):
            os.makedirs(vtk_dir)
        vtkFile = open(vtkFile1 + str(fr) + '.vtk', 'w')

        vtkFile.write('# vtk DataFile Version 2.0\n')
        vtkFile.write('ABAQUS Field Data\n')
        vtkFile.write('ASCII\n')
        vtkFile.write('DATASET UNSTRUCTURED_GRID\n')

        # Write points
        vtkFile.write('POINTS %i float\n' % numNodes)
        for nd in range(numNodes):
            x = InitialCoords[nd][0]
            y = InitialCoords[nd][1]
            z = InitialCoords[nd][2]
            if fr != 0 and fr in DISPLACEMENT:
                x += DISPLACEMENT[fr][nd].data[0]
                y += DISPLACEMENT[fr][nd].data[1]
                z += DISPLACEMENT[fr][nd].data[2]
            vtkFile.write('%f %f %f\n' % (x, y, z))

        # Write elements
        totalCellInts = sum([len(c)+1 for c in elementConnectivity])
        vtkFile.write('CELLS %i %i\n' % (numElements, totalCellInts))
        for conn in elementConnectivity:
            vtkFile.write('%i\t' % len(conn) + '\t'.join(str(n) for n in conn) + '\n')

        # Cell types
        vtkFile.write('CELL_TYPES %i\n' % numElements)
        for el in range(numElements):
            vtkFile.write('24\n')  # C3D10

        # -------------------------------
        # POINT_DATA (nodal)
        vtkFile.write('POINT_DATA %i\n' % numNodes)
        for fld in fieldDict.get(fr, {}):
            fieldValues = frames[fr].fieldOutputs[fld].getSubset(region=myInstance).values
            if is_nodal_field(fieldValues):
                # Determine if vector
                if isinstance(fieldValues[0].data, float):
                    # Scalar
                    vtkFile.write('SCALARS %s float 1\n' % fld)
                    vtkFile.write('LOOKUP_TABLE default\n')
                    for val in fieldDict[fr][fld]:
                        vtkFile.write('%f\n' % val)
                else:
                    # Vector
                    vtkFile.write('VECTORS %s float\n' % fld)
                    for val in fieldDict[fr][fld]:
                        vtkFile.write('%f %f %f\n' % (val[0], val[1], val[2]))

        # -------------------------------
        # CELL_DATA (element)
        vtkFile.write('CELL_DATA %i\n' % numElements)
        for fld in fieldDict.get(fr, {}):
            fieldValues = frames[fr].fieldOutputs[fld].getSubset(region=myInstance).values
            if not is_nodal_field(fieldValues):
                val0 = fieldDict[fr][fld][0]
                if isinstance(val0, float):
                    # Scalar
                    vtkFile.write('SCALARS %s float 1\n' % fld)
                    vtkFile.write('LOOKUP_TABLE default\n')
                    for val in fieldDict[fr][fld]:
                        vtkFile.write('%f\n' % val)
                elif len(val0) == 3:
                    # Vector
                    vtkFile.write('VECTORS %s float\n' % fld)
                    for val in fieldDict[fr][fld]:
                        vtkFile.write('%f %f %f\n' % (val[0], val[1], val[2]))
                elif len(val0) == 6:
                    # Symmetric tensor: Sxx,Syy,Szz,Sxy,Syz,Sxz
                    vtkFile.write('TENSORS %s float\n' % fld)
                    for val in fieldDict[fr][fld]:
                        Sxx,Syy,Szz,Sxy,Syz,Sxz = val
                        vtkFile.write('%f %f %f\n' % (Sxx,Sxy,Sxz))
                        vtkFile.write('%f %f %f\n' % (Sxy,Syy,Syz))
                        vtkFile.write('%f %f %f\n' % (Sxz,Syz,Szz))

        vtkFile.close()

    myOdb.close()

    '''
    abaqus viewer noGUI=odb_vtk_quad_tet_coarse_vb_nonlocal.py
    '''