import numpy as np
import pandas as pd
from pathlib import Path
import os

# Ring indices for DKP series
ring_indices = {
    'cGlyGly': [0, 1, 2, 3, 4, 5],  
    'cAlaAla': [4, 3, 2, 1, 0, 5],
    'cHisHis': [0, 1, 2, 3, 4, 5],
    'cHisHis2+': [0, 1, 2, 3, 4, 5],
    'cPhgPhg': [2, 1, 0, 5, 4, 3],
    'cLeuLeu': [0, 4, 5, 3, 2, 1],
    'cValVal': [4, 3, 2, 1, 0, 5],
    'cTrpTrp': [0, 1, 2, 3, 4, 5],
    'cPhePhe': [0, 5, 4, 3, 2, 1]
}

def parse_xyz_to_df(file_path):
    """
    Parse an XYZ file and return the number of atoms and a DataFrame with atomic coordinates.
    :param file_path: str, the path to the XYZ file to be parsed.
    :return: tuple, first element is an int representing the number of atoms,
             and the second element is a pandas DataFrame with columns 'Atom', 'X', 'Y', 'Z'.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()[2:]  # Skip the first two lines. 'lines' is a list of str.

        data = []  # This will be a list of lists, where each inner list represents an atom and its coordinates.

        for line in lines:
            parts = line.split()  # Splitting each line into parts. 'parts' is a list of str.
            atom = parts[0]  # 'atom' is a str representing the atomic symbol.
            x, y, z = map(float, parts[1:4])  # Converting strings to floats for the coordinates.
            data.append([atom, x, y, z])  # Appending the atom and its coordinates to the 'data' list.
    
    # Creating a DataFrame from the 'data' list. Each row represents an atom and its coordinates.
    df = pd.DataFrame(data, columns=['Atom', 'X', 'Y', 'Z'])
    num_atoms = len(df)
    return num_atoms, df

def get_ring_coord(df, indices):
    """
    Selects the atoms that compose the ring based on provided indices.
    """
    return df.iloc[indices]

def translate(coordinates):
    """
    Translate the ring coordinates to the origin center
    Input: coordinates (numpy array)
    Return: new_coordinates (numpy array)
    """
    return coordinates - coordinates.mean(axis=0)

def get_mean_plane(coordinates):
    """
    Compute the mean plane
    Input: coordinates (numpy array)
    Return: R1, R2 (numpy arrays)
    """
    N = coordinates.shape[0]  # Ring size
    R1 = np.dot(np.sin(2 * np.pi * np.arange(0, N) / N), coordinates)
    R2 = np.dot(np.cos(2 * np.pi * np.arange(0, N) / N), coordinates)
    return R1, R2

def get_normal(coordinates):
    """
    Compute normal to the mean plane
    Input: coordinates (numpy array)
    Output: unit_normal (numpy array)
    """
    R1, R2 = get_mean_plane(coordinates)
    cross_product = np.cross(R1, R2)
    return cross_product / np.linalg.norm(cross_product)

def fix_zero(x):
    """
    Fix the array at zero if the values are close enough to zero
    """
    return np.array([0.0]) if np.allclose(0, x, rtol=1e-06, atol=1e-08) else x

def displacement(coordinates):
    """
    Compute the displacement (z)
    Input: coordinates (numpy array)
    Output: Z (numpy array)
    """
    n = get_normal(coordinates)
    return np.dot(coordinates, n)

def get_ring_pucker_coords(coordinates):
    """
    Compute Ring Pucker Parameters
    Input: coordinates (numpy array)
    Return: puckering amplitude (q_i >= 0 for all i), angle in radians and degrees
    """
    N = coordinates.shape[0]  # Number of atoms in the ring
    z = displacement(coordinates)
    if 4 < N <= 20:
        if N % 2 == 0:
            m = range(2, int(N / 2))
            cos_component = [np.dot(z, np.cos(2 * np.pi * k * np.arange(0, N) / N)) for k in m]
            sin_component = [np.dot(z, np.sin(2 * np.pi * k * np.arange(0, N) / N)) for k in m]
            qcos = np.sqrt(2 / N) * np.array(cos_component)
            qsin = -np.sqrt(2 / N) * np.array(sin_component)
            q = np.sqrt(qsin**2 + qcos**2)
            amplitude = np.append(q, (1 / np.sqrt(N)) * np.dot(z, np.cos(np.arange(0, N) * np.pi)).sum())
            angle_rad = np.arctan2(qsin, qcos)
            angle_rad = np.where(angle_rad < 0, angle_rad + 2 * np.pi, angle_rad)
            angle_deg = np.degrees(angle_rad)
        else:
            m = range(2, int((N - 1) / 2) + 1)
            cos_component = [np.dot(z, np.cos(2 * np.pi * k * np.arange(0, N) / N)) for k in m]
            sin_component = [np.dot(z, np.sin(2 * np.pi * k * np.arange(0, N) / N)) for k in m]
            qcos = fix_zero(np.sqrt(2 / N) * np.array(cos_component))
            qsin = fix_zero(-np.sqrt(2 / N) * np.array(sin_component))
            amplitude = np.sqrt(qsin**2 + qcos**2)
            angle_rad = np.arctan2(qsin, qcos)
            angle_rad = np.where(angle_rad < 0, angle_rad + 2 * np.pi, angle_rad)
            angle_deg = np.degrees(angle_rad)
    else:
        print("Ring Size is too big or too small")
    return amplitude, angle_deg, angle_rad

def get_total_puckering_amplitude(amplitude):
    """
    Compute Total Puckering Amplitude
    Input: amplitude (numpy array)
    Output: Q (positive float)
    """
    return np.sqrt(np.square(amplitude).sum())

def get_spherical_polar_set_n6(amplitude, angle_deg):
    """
    Compute Spherical Polar Parameters (Q, θ, φ)
    Input: amplitude (numpy array), angle_deg (numpy array)
    Output: Q, theta_deg, phi_deg (floats)
    """
    Q = np.sqrt(np.square(amplitude).sum())
    q2 = amplitude[0]
    q3 = amplitude[1]
    
    theta_deg = np.degrees(np.arctan2(q2, q3))
    phi_deg = angle_deg[0]
    
    return Q, theta_deg, phi_deg

def conformation_haversine(amplitude, theta_deg, phi_deg):
    """
    Classifies the ring conformation based on θ and φ values using haversine.
    Returns the conformation as a string.
    """
    theta_deg = theta_deg % 180  # Normalize theta to [0, 180]
    phi_deg = phi_deg % 360      # Normalize phi to [0, 360]
    
    conformations = {
        'Chair (C)': {'theta': [0, 180]},
        'Boat (B)': {'theta': [90], 'phi': [0, 60, 120, 180, 240, 300]},
        'Twist-Boat (TB)': {'theta': [90], 'phi': [30, 90, 150, 210, 270, 330]},
        'Half-Chair (HC)': {'theta': [30, 150]},
        'Half-Boat (HB)': {'theta': [60, 120]}
    }
    
    min_distance = float('inf')
    assigned_conformation = None
    
    lat1 = np.deg2rad(90 - theta_deg)
    lon1 = np.deg2rad(phi_deg)
    point1 = np.array([[lat1, lon1]])
    
    for conformation, angles in conformations.items():
        theta_refs = np.deg2rad(angles.get('theta', [theta_deg]))
        phi_refs = np.deg2rad(angles.get('phi', [phi_deg]))
        
        for theta_ref in theta_refs:
            for phi_ref in phi_refs:
                lat2 = np.deg2rad(90 - np.degrees(theta_ref))
                lon2 = np.deg2rad(np.degrees(phi_ref))
                
                point2 = np.array([[lat2, lon2]])
                distance = np.linalg.norm(point1 - point2)
                
                if distance < min_distance:
                    min_distance = distance
                    assigned_conformation = conformation
                
    return assigned_conformation

def get_pucker_values(df, system_name):
    """
    Computes puckering values for a given system.
    :param df: DataFrame with atomic coordinates.
    :param system_name: str, name of the system being analyzed (e.g., 'cAlaAla').
    :return: puckering values and assigned conformation.
    """
    # Identify the ring atoms using the pre-defined ring indices for the system
    ring_atoms = get_ring_coord(df, ring_indices[system_name])
    
    # Translate ring atoms to center
    coordinates = translate(ring_atoms[['X', 'Y', 'Z']].values)
    
    # Get puckering amplitude, angles, and conformation
    amplitude, angle_deg, angle_rad = get_ring_pucker_coords(coordinates)
    Q, theta_deg, phi_deg = get_spherical_polar_set_n6(amplitude, angle_deg)
    
    conformation = conformation_haversine(amplitude, theta_deg, phi_deg)
    
    return {
        'amplitude': amplitude,
        'angle_deg': angle_deg,
        'Q': Q,
        'theta_deg': theta_deg,
        'phi_deg': phi_deg,
        'conformation': conformation
    }

def analyze_system(system_folder):
    results = []

    for chirality in ['SS', 'SR']:
        chirality_folder = Path(system_folder) / chirality
        print(f"Looking in folder: {chirality_folder}")

        if chirality_folder.exists():
            # Now recursively search for .xyz files in subdirectories
            xyz_files = list(chirality_folder.rglob('*.xyz'))
            print(f"  Found {len(xyz_files)} .xyz files")

            for xyz_file in xyz_files:
                system_name = xyz_file.parent.name  # Get system name from subfolder name
                file_stem = xyz_file.stem

                if system_name not in ring_indices:
                    print(f"Skipping {xyz_file.name}: system name '{system_name}' not in ring_indices")
                    continue

                try:
                    num_atoms, df = parse_xyz_to_df(xyz_file)
                    puckering_values = get_pucker_values(df, system_name)
                except Exception as e:
                    print(f"Error processing {xyz_file.name}: {e}")
                    continue

                result = {
                    'system_name': system_name,
                    'chirality': chirality,
                    'puckering_values': puckering_values
                }
                results.append(result)
        else:
            print(f"Folder does not exist: {chirality_folder}")
    
    return results

