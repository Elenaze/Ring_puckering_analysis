import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as Gridspec
import pandas as pd
import json
import os

def plot_puckering_distribution(json_path: str, save_path: str = None):
    """
    Plots puckering amplitude histogram and polar conformation plot
    from a JSON file containing puckering data.

    Parameters:
        json_path (str): Path to JSON file with puckering data
        save_path (str, optional): If provided, saves the figure to this path
    """

    # -----------------------------
    # 1. Load and parse JSON
    # -----------------------------
    if not os.path.isfile(json_path):
        print(f"JSON file not found: {json_path}")
        return

    with open(json_path, 'r') as f:
        raw_data = json.load(f)

    info = {
        key: pd.DataFrame([val])
        for key, val in raw_data.items()
    }

    # -----------------------------
    # 2. Plotting Setup
    # -----------------------------
    plt.rcParams['font.size'] = 16
    plt.rcParams['figure.dpi'] = 300

    fig = plt.figure(figsize=(8.5, 3.75))
    gs = Gridspec.GridSpec(1, 2, figure=fig, wspace=0.3)

    ax = [fig.add_subplot(gs[0, 0]),
          fig.add_subplot(gs[0, 1], projection='polar')]

    bin_width = 0.01
    bins = np.arange(0, 0.6, bin_width)

    capital_phi = "\u03A6"     # Φ
    capital_theta = "\u0398"   # Θ

    conformations = {
        'Chair (C)': {'theta': [0, 180]},
        'Boat (B)': {'theta': [90], 'phi': [0, 60, 120, 180, 240, 300]},
        'Twist-Boat (TB)': {'theta': [90], 'phi': [30, 90, 150, 210, 270, 330]},
        'Half-Chair (HC)': {'theta': [30, 150]},
        'Half-Boat (HB)': {'theta': [60, 120]},
    }

    marker_style = {
        'Chair (C)': {'marker': 'o', 'color': 'red'},
        'Boat (B)': {'marker': 's', 'color': 'green'},
        'Twist-Boat (TB)': {'marker': '^', 'color': 'orange'},
        'Half-Chair (HC)': {'marker': 'D', 'color': 'purple'},
        'Half-Boat (HB)': {'marker': 'v', 'color': 'blue'},
    }

    # -----------------------------
    # 3. Histogram and Polar Plot
    # -----------------------------
    for key, df_polar in info.items():
        ax[0].hist(df_polar['Amplitude'], bins=bins, alpha=0.7)

        theta_deg = df_polar['theta'] % 180
        phi_deg = df_polar['phi'] % 360
        phi_rad = np.deg2rad(phi_deg)

        ax[1].scatter(phi_rad, theta_deg, alpha=0.8, marker='.', s=250, edgecolors='black')

    ax[0].set_ylabel('Counts', color='black')
    ax[0].set_xlabel(r'Q / $\AA$', color='black')
    ax[0].set_xticks([0.0, 0.2, 0.4, 0.6])
    ax[0].set_yticks([0, 2.5, 5, 7.5, 10])
    ax[0].tick_params(axis='both', colors='black', width=1.)

    ax[1].set_ylim(0, 180)

    for conformation, angles in conformations.items():
        theta_values = angles.get('theta', [])
        phi_values = angles.get('phi', [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
        style = marker_style[conformation]

        for theta in theta_values:
            for phi in phi_values:
                theta_rad = np.deg2rad(phi)
                r = theta
                ax[1].plot(theta_rad, r, marker=style['marker'], markersize=5, alpha=0.1,
                           color=style['color'], linestyle='None', label=conformation)
                
    #Legend if needed uncomment the following lines
    #--------------------------------------------
    #handles, labels = ax[1].get_legend_handles_labels()
    #by_label = dict(zip(labels, handles))
    #ax[1].legend(by_label.values(), by_label.keys(), loc='upper right', bbox_to_anchor=(1.15, 1.05))

    ax[1].set_xticks(np.deg2rad([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]))
    ax[1].set_xticklabels(['0°', '30°', '60°', '90°', '120°', '150°', '180°', '210°', '240°', '270°', '300°', '330°'], fontsize=16, color='black')

    ax[1].set_yticklabels([])
    ax[1].set_xlabel(f"{capital_phi} / °", color='black')
    ax[1].set_ylabel(f"{capital_theta} / °", color='black', labelpad=30)
    ax[1].yaxis.set_label_position("right")
    ax[1].tick_params(pad=8)

    if save_path:
        fig.savefig(save_path, bbox_inches='tight', transparent=True, dpi=300)
        print(f"Plot saved to: {save_path}")

    plt.show()
