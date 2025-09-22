# Dipeptide Ring Puckering Analysis

This analysis tool is specifically designed for the research article on cyclic dipeptide (diketopiperazine, DKP) conformational analysis. It provides specialized tools for analyzing the puckering parameters of DKP rings from the computational study.

**Note**: For general-purpose chirality analysis tools, see the `chirality_descriptors/` directory which contains reusable modules for various molecular systems.

## Overview

This article-specific analysis characterizes the three-dimensional shape of cyclic dipeptides by computing puckering parameters including:
- **Total puckering amplitude (Q)**: Overall degree of ring puckering
- **Spherical polar coordinates (θ, φ)**: Angular parameters describing puckering geometry
- **Conformation classification**: Assignment to standard ring conformations (Chair, Boat, Twist-Boat, etc.)

## Features

- **XYZ file parsing**: Reads molecular coordinates from standard XYZ format files
- **Automated ring detection**: Uses predefined atom indices for different DKP systems
- **Puckering parameter calculation**: Computes Cremer-Pople puckering parameters
- **Conformation assignment**: Classifies rings into standard conformational categories
- **Data visualization**: Generates histograms and polar plots of puckering distributions
- **Batch processing**: Analyzes multiple structures across different chiralities

## Supported Systems

The following cyclic dipeptide systems are supported with predefined ring atom indices:

- `cGlyGly` - Cyclo(Gly-Gly)
- `cAlaAla` - Cyclo(Ala-Ala) 
- `cHisHis` - Cyclo(His-His)
- `cHisHis2+` - Cyclo(His-His) doubly protonated
- `cPhgPhg` - Cyclo(Phg-Phg)
- `cLeuLeu` - Cyclo(Leu-Leu)
- `cValVal` - Cyclo(Val-Val)
- `cTrpTrp` - Cyclo(Trp-Trp)
- `cPhePhe` - Cyclo(Phe-Phe)

## Related Tools

This analysis is part of a larger study that also uses general-purpose chirality descriptors:
- **Alpha determination**: PDF fitting analysis (see `chirality_descriptors/alpha_determination/`)
- **mo_RMSD**: Modified RMSD chirality descriptor (see `chirality_descriptors/mo_RMSD/`)

## Installation

Ensure you have the required Python packages installed:

```bash
pip install numpy pandas matplotlib pathlib
```

## Usage

### Basic Usage

Run the analysis from the command line:

```bash
python main.py <path_to_data_directory>
```

Example:
```bash
python main.py ../data
```

### Expected Directory Structure

Your data directory should be organized as follows:

```
data/
├── SS/
│   ├── cAlaAla/
│   │   └── structure.xyz
│   ├── cGlyGly/
│   │   └── structure.xyz
│   └── ...
└── SR/
    ├── cAlaAla/
    │   └── structure.xyz
    ├── cGlyGly/
    │   └── structure.xyz
    └── ...
```

Where:
- `SS/` and `SR/` represent different chiralities
- Each system folder contains XYZ coordinate files
- XYZ files should follow standard format (atom count, comment line, then atomic coordinates)

### Output

The analysis generates:

1. **JSON data file** (`output/puckering_data.json`): Contains all computed puckering parameters
2. **Summary plot** (`output/puckering_summary.pdf`): Visualization showing:
   - Histogram of puckering amplitudes
   - Polar plot of conformational distribution

## Core Functions

### `ring_analysis.py`

#### Key Functions:

- `parse_xyz_to_df(file_path)`: Parses XYZ files into pandas DataFrames
- `get_pucker_values(df, system_name)`: Computes all puckering parameters for a structure
- `get_ring_pucker_coords(coordinates)`: Calculates Cremer-Pople puckering coordinates
- `conformation_haversine(amplitude, theta_deg, phi_deg)`: Assigns conformational labels
- `analyze_system(system_folder)`: Batch processes all structures in a directory

#### Ring Conformations:

The script classifies rings into the following conformations based on θ and φ angles:

- **Chair (C)**: θ ≈ 0° or 180°
- **Boat (B)**: θ ≈ 90°, φ = 0°, 60°, 120°, 180°, 240°, 300°
- **Twist-Boat (TB)**: θ ≈ 90°, φ = 30°, 90°, 150°, 210°, 270°, 330°
- **Half-Chair (HC)**: θ ≈ 30° or 150°
- **Half-Boat (HB)**: θ ≈ 60° or 120°

![Polar Plot Conformations](polar_plot_conf.pdf)

*Figure: Polar coordinate representation showing the angular regions corresponding to different ring conformations. The radial axis represents θ (0° to 180°) and the angular axis represents φ (0° to 360°). Different markers indicate the canonical positions for each conformation type.*

### `plotting.py`

Generates plots including:
- Puckering amplitude histograms
- Polar coordinate conformational maps
- Reference conformational markers

## Theory

The puckering analysis is based on the Cremer-Pople method for characterizing ring conformations:

1. **Coordinate transformation**: Ring atoms are translated to center and projected onto mean plane
2. **Displacement calculation**: Out-of-plane displacements (z-coordinates) are computed
3. **Puckering coordinates**: Amplitude and phase parameters are derived from Fourier analysis
4. **Spherical coordinates**: Total amplitude Q and angular parameters θ, φ are calculated

For 6-membered rings:
- **Q**: Total puckering amplitude = √(q₂² + q₃²)
- **θ**: Polar angle describing chair vs. non-chair character
- **φ**: Azimuthal angle describing specific non-chair conformations

## File Formats

### Input XYZ Format:
```
6
Comment line
C    x1    y1    z1
N    x2    y2    z2
...
```

### Output JSON Format:
```json
{
    "cAlaAla_SS": {
        "Amplitude": 0.123,
        "theta": 45.6,
        "phi": 78.9,
        "conformation": "Half-Chair (HC)"
    }
}
```

## Error Handling

The script includes robust error handling for:
- Missing or invalid XYZ files
- Unsupported system names
- Invalid coordinate data
- File I/O errors

## References

- Cremer, D. & Pople, J. A. (1975). General definition of ring puckering coordinates. *J. Am. Chem. Soc.*, 97, 1354-1358.
- Boeyens, J. C. A. (1978). The conformation of six-membered rings. *J. Cryst. Mol. Struct.*, 8, 317-320.

## License

This software is provided as-is for research purposes.