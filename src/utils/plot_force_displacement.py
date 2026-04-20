'''
This script can be used to create force-displacement plots from FE simulations using the non-local model.
'''

#########################################################################################################

## Modules
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

## User Input
# Model data
folder_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/'
model_name = 'C3D10_cube_NL.odb'   
step_name = 'Step-1'      
out_name = 'RF_history.csv'
top_nodes = [5,6,7,8,13,14,15,16,22]

# Other input
utils_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/src/utils/'

#**************************************************************************************************************

## Extract data from odb
odb_path = folder_path + model_name
out_csv = folder_path + out_name

cmd = ["abaqus", "python", utils_path + "extract_RF_from_ODB.py", odb_path, step_name, out_csv]
result = subprocess.run(cmd)

if result.returncode != 0:
        print("Abaqus python data extraction failed:")
        print(result.stderr)
else:
        print("Extraction finished successfully")

## Import results
df_res = pd.read_csv(out_csv)

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
plotname = folder_path + model_name[:-3] + "_fd.png"
plt.figure()
plt.plot(res['U3'], res['RF3_sum'])
plt.xlabel('displacement')
plt.ylabel('force')
plt.savefig(plotname, dpi=300 )
plt.show()