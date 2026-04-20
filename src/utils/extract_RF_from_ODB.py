from odbAccess import openOdb
import csv

odb_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/C3D10_cube_NL.odb'      # <-- change
step_name = 'Step-1'      # <-- change
out_csv = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/RF_history.csv'

odb = openOdb(odb_path)
step = odb.steps[step_name]

f = open(out_csv, 'wb')   # Python 2 requires binary mode
writer = csv.writer(f)
writer.writerow(['Increment', 'FrameValue', 'Node', 'RF1', 'RF2', 'RF3'])

for frame in step.frames:
    inc = frame.incrementNumber
    time = frame.frameValue

    if 'RF' not in frame.fieldOutputs:
        continue

    rf_field = frame.fieldOutputs['RF']

    for v in rf_field.values:
        node = v.nodeLabel
        rf = v.data

        rf1 = rf[0]
        rf2 = rf[1] if len(rf) > 1 else 0.0
        rf3 = rf[2] if len(rf) > 2 else 0.0

        writer.writerow([inc, time, node, rf1, rf2, rf3])

f.close()
odb.close()