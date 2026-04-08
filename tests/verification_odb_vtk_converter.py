import numpy as np
import meshio
import glob
import os
import re

# -------- USER INPUT --------
dat_file = "/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/single_C3D10_Abq/C3D10_single_Abq.dat"
vtk_folder = "/home/gabriela/Documents/04_Projects/2026_NonLocal_Damage_Model/02_Code/non-local-damage-model/tests/single_C3D10_Abq/vtk_C3D10_single_Abq/"
vtk_pattern = "C3D10_single_Abq*.vtk"
tol = 1e-6
# ----------------------------

# ---------------- DAT PARSING ----------------
import numpy as np
import re
from collections import defaultdict


def parse_dat_increment_aligned(dat_file):
    """Parse Abaqus DAT file and align outputs per increment."""

    with open(dat_file, "r") as f:
        lines = f.readlines()

    increments = []
    current_inc = {}
    increment_active = False

    # --------- YOUR ORIGINAL TABLE PARSER (unchanged) ----------
    def parse_table_block(start_idx):
        node_ids = []
        data = []
        i = start_idx
        while i < len(lines):
            row = lines[i].strip()
            if row == "" or "MAXIMUM" in row or "MINIMUM" in row:
                break
            parts = row.split()
            try:
                node_ids.append(int(parts[0]))
                if len(parts) > 4:
                    values = [float(x) for x in parts[1:]]
                else:
                    values = [float(x) for x in parts]
                data.append(values[1:])
            except:
                break
            i += 1
        return i, np.array(node_ids), np.array(data)

    # --------- NEW SDV BLOCK PARSER ----------
    def parse_sdv_block(start_idx, sdv_store):
        """
        Parse one SDV table and append its columns to sdv_store[(elem,pt)]
        """
        header_line = lines[start_idx + 3]
        sdv_ids = [int(x.replace("SDV", "")) for x in re.findall(r"SDV\d+", header_line)]

        i = start_idx + 6  # where data rows start

        while i < len(lines):
            row = lines[i].strip()

            if row == "" or "MAXIMUM" in row:
                break

            parts = row.split()

            try:
                elem = int(parts[0])
                pt = int(parts[1])
                values = [float(x) for x in parts[2:2+len(sdv_ids)]]
            except:
                break

            sdv_store[(elem, pt)].extend(values)
            i += 1

        return i

    # ------------------------------------------------------------

    i = 0
    sdv_store = None

    while i < len(lines):
        line = lines[i]

        # --------- Detect new increment ----------
        if "STEP" in line and "INCREMENT" in line:
            if increment_active:
                # finalize SDVs for previous increment
                if sdv_store:
                    keys = sorted(sdv_store.keys())
                    sdv_array = np.array([sdv_store[k] for k in keys])
                    current_inc["SDV"] = sdv_array

                increments.append(current_inc)

            current_inc = {"U": None, "RF": None, "SDV": None, "LE": None, "S": None}
            sdv_store = defaultdict(list)
            increment_active = True

        # --------- U block ----------
        if ("THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ALLNODES" in line
                and "U1" in lines[i+2]):
            i, node_ids, data = parse_table_block(i+5)
            current_inc["U"] = data
            continue

        # --------- RF block ----------
        if ("THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ALLNODES" in line
                and "RF1" in lines[i+2]):
            i, node_ids, data = parse_table_block(i+5)
            current_inc["RF"] = data
            continue

        # --------- SDV block ----------
        if ("THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS" in line
                and "SDV" in lines[i+3]):
            i = parse_sdv_block(i, sdv_store)
            continue

        # --------- LE block ----------
        if ("THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS" in line
                and re.search(r"E1\d", lines[i+3])):
            i, node_ids, data = parse_table_block(i+6)
            current_inc["LE"] = data
            continue

        # --------- S block ----------
        if ("THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS" in line
                and re.search(r"S1\d", lines[i+3])):
            i, node_ids, data = parse_table_block(i+6)
            current_inc["S"] = data
            continue

        i += 1

    # --------- Finalize last increment ----------
    if increment_active:
        if sdv_store:
            keys = sorted(sdv_store.keys())
            sdv_array = np.array([sdv_store[k] for k in keys])
            current_inc["SDV"] = sdv_array

        increments.append(current_inc)

    return increments

# ---------------- VTK PARSING ----------------
def read_vtk_file(vtk_file):
    """Read VTK file and separate point and cell data."""
    vtk = meshio.read(vtk_file)
    point_data = {
        "U": vtk.point_data.get("U"),
        "RF": vtk.point_data.get("RF")
    }
    cell_data = {
        "LE": vtk.cell_data_dict.get("LE"),
        "S": vtk.cell_data_dict.get("S"),
        "SDV": {k: vtk.cell_data_dict.get(k) for k in vtk.cell_data_dict if k.startswith("SDV")}
    }
    return {"points": vtk.points, "point_data": point_data, "cell_data": cell_data}

# ---------------- UTILITY FUNCTIONS ----------------
def reshape_dat_to_vtk(dat_arr, vtk_arr_shape):
    """
    Average DAT over integration points and convert Voigt 6 -> 3x3 if needed.
    VTK has one tensor per element.
    """
    if dat_arr is None or vtk_arr_shape is None:
        return None

    dat_arr = np.array(dat_arr)

    # Average over integration points
    dat_arr = np.mean(dat_arr, axis=0)  # (n_elem, n_vals)

    # Convert Voigt 6 -> 3x3 if VTK expects 3x3 tensor
    if dat_arr.shape[0] == 6 and len(vtk_arr_shape) == 3 and vtk_arr_shape[2] == 3:
        dat_arr = voigt6_to_3x3(dat_arr)

    return dat_arr

def voigt6_to_3x3(dat_voigt):
    """
    Convert an array in Voigt 6 notation to 3x3 symmetric matrices.
    dat_voigt: (n_elements, n_ip, 6)
    Returns: (n_elements, n_ip, 3, 3)
    """
    dat_voigt = np.array(dat_voigt)
    dat_3x3 = np.zeros((3, 3))
    dat_3x3[0, 0] = dat_voigt[0]  # S11
    dat_3x3[1, 1] = dat_voigt[1]  # S22
    dat_3x3[2, 2] = dat_voigt[2]  # S33
    dat_3x3[0, 1] = dat_voigt[3]  # S12
    dat_3x3[1, 0] = dat_voigt[3]  # symmetric
    dat_3x3[0, 2] = dat_voigt[4]  # S13
    dat_3x3[2, 0] = dat_voigt[4]  # symmetric
    dat_3x3[1, 2] = dat_voigt[5]  # S23
    dat_3x3[2, 1] = dat_voigt[5]  # symmetric
    return dat_3x3

def compare_arrays(dat_arr, vtk_arr, tol=1e-6, name="field"):
    """Compare two arrays and print summary statistics."""
    if vtk_arr is None or dat_arr is None:
        print(f"{name} missing, skipping")
        return
    if dat_arr.shape != vtk_arr.shape:
        print(f"{name}: shape mismatch DAT {dat_arr.shape} vs VTK {vtk_arr.shape}")
        return
    diff = vtk_arr - dat_arr
    norm_diff = np.linalg.norm(diff, axis=1) if diff.ndim == 2 else np.abs(diff)
    max_diff = np.max(norm_diff)
    mean_diff = np.mean(norm_diff)
    n_fail = np.sum(norm_diff > tol)
    print(f"{name}: max diff = {max_diff:.3e}, mean diff = {mean_diff:.3e}, exceed tol = {n_fail}")

# ---------------- MAIN SCRIPT ----------------
dat_increments = parse_dat_increment_aligned(dat_file)
vtk_files = sorted(glob.glob(os.path.join(vtk_folder, vtk_pattern)))

print(f"Found {len(dat_increments)} increments in DAT")
print(f"Found {len(vtk_files)} VTK files")

for i, vtk_file in enumerate(vtk_files):
    print(f"\n=== Increment {i+1} / {len(dat_increments)}: {vtk_file} ===")
    vtk = read_vtk_file(vtk_file)

    if i == 0:
        continue  # skip initial increment

    if i >= len(dat_increments):
        print("No corresponding DAT increment found")
        continue

    dat_inc = dat_increments[i]

    # Compare point data
    for field in ["U", "RF"]:
        vtk_arr = vtk["point_data"].get(field)
        compare_arrays(dat_inc[field], vtk_arr, tol, field)

    # Compare cell data (LE, S)
    for field in ["LE", "S"]:
        cell_dict = vtk["cell_data"].get(field)
        if cell_dict is None:
            vtk_arr = None
        else:
            vtk_arr = list(cell_dict.values())[0]
        dat_arr_expanded = reshape_dat_to_vtk(dat_inc[field], vtk_arr.shape)
        vtk_arr = vtk_arr[0]
        compare_arrays(dat_arr_expanded, vtk_arr, tol, field)

    # Compare SDV fields
    if dat_inc["SDV"] is not None and vtk["cell_data"]["SDV"]:
        for j, sdv_name in enumerate(sorted(vtk["cell_data"]["SDV"].keys())):
            vtk_sdv_dict = vtk["cell_data"]["SDV"][sdv_name]
            vtk_sdv = list(vtk_sdv_dict.values())[0]
            dat_sdv_expanded = reshape_dat_to_vtk(dat_inc["SDV"][:, j][:, np.newaxis], vtk_sdv.shape)
            compare_arrays(dat_sdv_expanded, vtk_sdv[0], tol, sdv_name)