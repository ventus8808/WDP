import geopandas as gpd
import pandas as pd
import yaml
from pathlib import Path

def create_adjacency_matrix(shapefile_path, output_path):
    """
    Creates a county adjacency matrix from a shapefile and saves it as a CSV.

    The matrix will have county GEOID as both the index and columns.
    A value of True indicates that two counties are adjacent (i.e., their boundaries touch).

    Args:
        shapefile_path (str or Path): Path to the county shapefile.
        output_path (str or Path): Path to save the output CSV file.
    """
    try:
        # Load the shapefile
        print(f"Loading shapefile from {shapefile_path}...")
        counties = gpd.read_file(shapefile_path)
        print(f"Shapefile loaded successfully with {len(counties)} counties.")

        # Ensure the required 'GEOID' column exists for county FIPS codes
        if 'GEOID' not in counties.columns:
            raise ValueError("Shapefile is missing the required 'GEOID' column.")

        # Use GEOID as the unique identifier for the matrix
        counties = counties.set_index('GEOID')

        # Create an empty DataFrame to store the adjacency matrix, initialized with False
        print("Initializing adjacency matrix...")
        adj_matrix = pd.DataFrame(False, index=counties.index, columns=counties.index)

        # Determine adjacency for each county
        print("Calculating adjacencies... This may take a few minutes.")
        for index, county in counties.iterrows():
            # Find all geometries that touch the current county's geometry
            neighbors = counties[counties.geometry.touches(county.geometry)].index.tolist()

            # Set the corresponding entry in the matrix to True for each neighbor
            if neighbors:
                adj_matrix.loc[index, neighbors] = True

        print("Adjacency matrix created successfully.")

        # Create edge list format for Bayesian models
        print("Creating edge list format...")
        edge_list = []
        for county_from in adj_matrix.index:
            for county_to in adj_matrix.columns:
                if adj_matrix.loc[county_from, county_to] == True:
                    edge_list.append({
                        'county_from': county_from,
                        'county_to': county_to,
                        'adjacency_weight': True
                    })

        edge_df = pd.DataFrame(edge_list)
        print(f"Created edge list with {len(edge_df)} adjacency relationships")

        # Ensure the output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save the matrix to a CSV file
        adj_matrix.to_csv(output_path)
        print(f"Adjacency matrix saved to {output_path}")

        # Save the edge list to a separate CSV file
        edge_list_path = output_path.parent / 'County_Adjacency_List.csv'
        edge_df.to_csv(edge_list_path, index=False)
        print(f"Adjacency edge list saved to {edge_list_path}")

    except FileNotFoundError:
        print(f"Error: Shapefile not found at {shapefile_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def main():
    """
    Main function to orchestrate the script execution.
    It loads configuration, defines paths, and calls the matrix creation function.
    """
    # Define the project root relative to this script's location
    # WDP/Code/Clean/County_Adjacency.py -> WDP/
    project_root = Path(__file__).resolve().parents[2]
    
    # Load the project configuration file
    config_path = project_root / 'config.yaml'
    
    try:
        print(f"Loading configuration from {config_path}...")
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Get the relative paths from the configuration file
        shapefile_rel_path = config['data_sources']['tiger']['shapefile']
        socioeconomic_dir_rel_path = config['data_sources']['socioeconomic']['saipe']['processed']

        # Construct the full, absolute paths
        shapefile_path = project_root / shapefile_rel_path
        output_dir = project_root / socioeconomic_dir_rel_path
        output_file_path = output_dir / 'County_Adjacency_Matrix.csv'

        # Generate and save the adjacency matrix
        create_adjacency_matrix(shapefile_path, output_file_path)

    except FileNotFoundError:
        print(f"Error: Configuration file not found at {config_path}")
    except KeyError as e:
        print(f"Error: Missing expected key in config.yaml: {e}")
    except Exception as e:
        print(f"An unexpected error occurred in main(): {e}")

if __name__ == "__main__":
    main()
