'''
This script can be used to create force-displacement plots from FE simulations using the non-local model.
'''

#########################################################################################################

## Modules
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

#**************************************************************************************************************

## Extract data from odb
def extract_fd_from_odb(folder_path,model_name,out_name):
    odb_path = folder_path + model_name
    out_csv = folder_path + out_name

    cmd = ["abaqus", "python", utils_path + "extract_RF_from_ODB.py", odb_path, step_name, out_csv]
    result = subprocess.run(cmd)

    if result.returncode != 0:
            print("Abaqus python data extraction failed:")
            print(result.stderr)
    else:
            print("Extraction finished successfully")

# Filter force displacement data
def get_fd_curve(path_results, top_nodes):
    ## Import results
    df_res = pd.read_csv(path_results)

    # Ensure numeric colums are numeric
    cols = ['Increment','RF3','U3','Node']
    df_res[cols] = df_res[cols].apply(pd.to_numeric)

    ## Calculate reaction force
    results = []

    for inc, g in df_res.groupby('Increment'):

        gsel = g[g['Node'].isin(top_nodes)]

        Fz = gsel['RF3'].sum()
        Uz = gsel['U3'].mean()

        results.append([inc, Uz, Fz])

    res = pd.DataFrame(results, columns=['Increment','U3','RF3_sum'])

    return res['U3'].values, res['RF3_sum'].values

# Make force-displacement plot
def fd_plot(plotname,f,d):
    plt.figure()
    plt.plot(d, f)
    plt.xlabel('displacement')
    plt.ylabel('force')
    plt.savefig(plotname, dpi=300 )
    plt.show()  

# Force-displacement plot overlay function
def fd_plot_overlay(plotname, curves, labels):
    plt.figure()

    for (d, f), label in zip(curves, labels):
        plt.plot(d, f, label=label)

    plt.xlabel('displacement')
    plt.ylabel('force')
    plt.legend()
    plt.grid(True)
    plt.savefig(plotname, dpi=300)
    plt.show()

if __name__ == "__main__":
    
    ## User Input
    # Model data
    path1 = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/'
    path2 = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_L_comp/'
    path3 = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_Abq_comp/'

    model_name1 = 'C3D10_cube_NL.odb'
    model_name2 = 'C3D10_cube_NL.odb'
    model_name3 = 'C3D10_cube_Abq.odb'

    model1_label = 'Nonlocal'
    model2_label = 'Local'
    model3_label = 'Abaqus'

    step_name = 'Step-1'      
    out_name = 'RF_history.csv'
    top_nodes = [5,6,7,8,13,14,15,16,22]

    # Other input
    utils_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/src/utils/'

    # Assemble input information
    folder_path = [path1, path2, path3]
    model_name = [model_name1, model_name2, model_name3]  
    labels = [model1_label, model2_label, model3_label]
    curves = []

    for fp, mn, label in zip(folder_path, model_name, labels):

        # Extract data from odb
        extract_fd_from_odb(fp, mn, out_name)

        csv_path = fp + out_name

        # Compute curve
        d, f = get_fd_curve(csv_path, top_nodes)
        curves.append((d, f))

        # Individual plot
        plotname = fp + mn[:-4] + "_fd.png"
        fd_plot(plotname, f, d)

    # --- Overlay plot
    plotname_overlay = folder_path[0] + "FD_overlay.png"
    fd_plot_overlay(plotname_overlay, curves, labels)

        