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

Supporting data:
- inhomCube_C3D10_L
- inhomCube_C3D10_NL_xl03
- inhomCube_C3D10_NL
- inhomCube_C3D10_NLxl15

The larger the length scale, the more homogeneous is the damage.
This also has an impact on the deformation on the object.

****************************************************************

VERIFICATION 4

****************************************************************

Objective:
When using the UEL with l_nl=0, i.e. a local damage formulation, 
the results should be the same as when using no UEL and just 
having a UMAT without any non-local parameters.

Supporting data:
- homCube_C3D10_L_comp
- homCube_C3D10_Abq_comp

The results in the dat file between the two simulations are 
identical. However, when printing dammage values from the 
UMAT, small differences (< 10^-5) can be observed, indicating 
that the results are not perfectly equaivalent. However, it 
should be noted that the UMAT had to be modified in order to
run it independently from the UEL. Furthermore all non-local
parameters were removed from the UMAT code.
