import os
import json
import pandas as pd
from pathlib import Path
import sys
from ring_analysis import analyze_system
from plotting import plot_puckering_distribution

# Get the directory containing main.py
SCRIPT_DIR = Path(__file__).parent.absolute()

def print_summary(results):
    """
    Print a summary of the puckering analysis for each system.
    """
    print("Puckering Analysis Summary:")
    for result in results:
        system_name = result['system_name']
        chirality = result['chirality']
        puckering_values = result['puckering_values']
        conformation = puckering_values['conformation']
        Q = puckering_values['Q']
        theta_deg = puckering_values['theta_deg']
        phi_deg = puckering_values['phi_deg']
        
        print(f"System: {system_name} ({chirality})")
        print(f"  Conformation: {conformation}")
        print(f"  Total Puckering Amplitude (Q): {Q:.3f}")
        print(f"  Theta (deg): {theta_deg:.2f}")
        print(f"  Phi (deg): {phi_deg:.2f}")
        print('-' * 40)


def save_results(results, output_dir):
    """Save results to JSON file in the specified output directory"""
    output_file = output_dir / "puckering_data.json"
    
    export_dict = {}
    for result in results:
        system = result['system_name']
        chirality = result['chirality']
        key = f"{system}_{chirality}"
        
        puck = result['puckering_values']
        export_dict[key] = {
            'Amplitude': puck['Q'],
            'theta': puck['theta_deg'],
            'phi': puck['phi_deg'],
            'conformation': puck['conformation']
        }

    with open(output_file, 'w') as f:
        json.dump(export_dict, f, indent=4)

    print(f"Full puckering data saved to: {output_file}")

def main() -> bool:
    """Main function to run puckering analysis"""
    if len(sys.argv) != 2:
        print("Usage: python main.py <data_directory>")
        sys.exit(1)
    
    # Setup paths
    data_path = SCRIPT_DIR / 'data'
    
    output_dir = SCRIPT_DIR / 'output'
    output_dir.mkdir(exist_ok=True)

    try:
        # Validate input directory
        if not data_path.is_dir():
            raise FileNotFoundError(f"Directory {data_path} does not exist.")
        
        # Run analysis
        results = analyze_system(str(data_path))
        if not results:
            raise ValueError("No results found in the analysis.")
        
        # Process and save results
        #print_summary(results)
        save_results(results, output_dir)
        
        # Generate plots
        plot_puckering_distribution(
            json_path=str(output_dir / "puckering_data.json"),
            save_path=str(output_dir / "puckering_summary.pdf")
        )
        
        return True

    except FileNotFoundError as e:
        print(f"File system error: {str(e)}")
        return False
    except ValueError as e:
        print(f"Value error: {str(e)}")
        return False
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        return False

if __name__ == "__main__":
    success = main()
    if not success:
        print("Analysis failed. Please check the error messages above.")