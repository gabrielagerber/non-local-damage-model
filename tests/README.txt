******************************************************************

VERIFICATION 1

******************************************************************

Objective:
In case of homogeneous stress and strain states, the damage is 
also homogeneous for all elements. In such a case, there is no 
difference between the local and the non-local damage model.

Supporting data:
- homCube_C3D10_L_comp
- homCube_C3D10_NL_comp

The results from these simulations are identical.

******************************************************************

VERIFICATION 2

******************************************************************

Objective:
For inhomogeneous stress and strain states, differences in the damage
variable between elements are expected. Hence, there should be a 
difference in the simulation results between local and non-local
simulations.

Supporting data:
- inhomCube_C3D10_L
- inhomCube_C3D10_NL

The damage values in the local and in the non-local simulation are
very different. The local formulation leads to a smoother damage 
distribution.

******************************************************************

VERIFICATION 3

******************************************************************

Objective:
With decreasing length-scale parameter XL, the impact of the non-
local damage model should become smaller and smaller.

Supporring data:
-
- 

ADD CONCLUSION HERE.
