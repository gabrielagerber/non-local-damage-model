import sys

def FD_from_odb(odb_path, step_name, out_csv):
    from odbAccess import openOdb
    import csv
    '''
    This function extracts force and displacement information for each node from .odf and stores
    it in a csv file.
    '''

    odb = openOdb(odb_path)
    step = odb.steps[step_name]

    f = open(out_csv, 'wb')   # Python 2 requires binary mode
    writer = csv.writer(f)

    writer.writerow([
        'Increment', 'FrameValue', 'Node',
        'RF1', 'RF2', 'RF3',
        'U1', 'U2', 'U3'
    ])

    for frame in step.frames:
        inc = frame.incrementNumber
        time = frame.frameValue

        # Skip frame if fields are missing
        if 'RF' not in frame.fieldOutputs or 'U' not in frame.fieldOutputs:
            continue

        rf_field = frame.fieldOutputs['RF']
        u_field  = frame.fieldOutputs['U']

        # Build fast lookup for displacement by node label
        u_dict = {}
        for v in u_field.values:
            u_dict[v.nodeLabel] = v.data

        for v in rf_field.values:
            node = v.nodeLabel
            rf = v.data

            # Reaction forces
            rf1 = rf[0]
            rf2 = rf[1] if len(rf) > 1 else 0.0
            rf3 = rf[2] if len(rf) > 2 else 0.0

            # Displacements (match by node)
            if node in u_dict:
                u = u_dict[node]
                u1 = u[0]
                u2 = u[1] if len(u) > 1 else 0.0
                u3 = u[2] if len(u) > 2 else 0.0
            else:
                u1 = u2 = u3 = 0.0

            writer.writerow([inc, time, node, rf1, rf2, rf3, u1, u2, u3])

    f.close()
    odb.close()

if __name__ == "__main__":
    # odb_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/C3D10_cube_NL.odb'   
    # step_name = 'Step-1'      
    # out_csv = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/RF_history.csv'
    
    odb_path  = sys.argv[1]
    step_name = sys.argv[2]
    out_csv   = sys.argv[3]

    FD_from_odb(odb_path, step_name, out_csv)