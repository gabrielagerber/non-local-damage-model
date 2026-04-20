'''
This script can be used to create force-displacement plots from FE simulations using the non-local model.
'''

#########################################################################################################

## Modules
from src.utils import extract_RF_from_ODB

## User Input
folder_path = '/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/homCube_C3D10_NL_comp/'
model_name = 'C3D10_cube_NL.odb'   
step_name = 'Step-1'      
out_name = 'RF_history.csv'

## Extract data from odb
odb_path = folder_path + model_name
out_csv = folder_path + out_name

## Import results

## Calculate reaction force

## Plot