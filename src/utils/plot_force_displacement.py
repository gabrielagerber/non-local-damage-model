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

def fd_plot(plotname,f,d):
    plt.figure()
    plt.plot(d, f)
    plt.xlabel('displacement')
    plt.ylabel('force')
    plt.savefig(plotname, dpi=300 )
    plt.show()
                 

## Create Force-displacement plot
def create_fd_plot(path_results, path, model_name):

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

    ## Plot
    plotname = path + model_name[:-4] + "_fd.png"
    f = res['RF3_sum']
    d = res['U3']
    fd_plot(plotname,f,d)
    

if __name__ == "__main__":
    
    ## User Input
    # Model data
    path1 = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/'
    path2 = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_L_comp/'
    path3 = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_Abq_comp/'

    model_name1 = 'C3D10_cube_NL.odb'
    model_name2 = 'C3D10_cube_NL.odb'
    model_name3 = 'C3D10_cube_Abq.odb'

    folder_path = [path1, path2, path3]
    model_name = [model_name1, model_name2, model_name3]   
    step_name = 'Step-1'      
    out_name = 'RF_history.csv'
    top_nodes = [5,6,7,8,13,14,15,16,22]

    # Other input
    utils_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/src/utils/'

    # Extract results
    for n, path in enumerate(folder_path):
        extract_fd_from_odb(path,model_name[n],out_name)

        # Create force-displacement plot
        path_results = path + out_name
        create_fd_plot(path_results, path, model_name[n])

        