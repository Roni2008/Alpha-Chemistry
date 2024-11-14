
import zipfile
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os
import csv
from mpl_toolkits.mplot3d import Axes3D  
from matplotlib.widgets import Button
import tkinter as tk
from tkinter import simpledialog, messagebox
import Alpha11
from Alpha11 import *
import os
from tkinter import Tk, filedialog
import networkx as nx

def get_connected_nodes(bonds_df, from_node, to_node):
    # Create a graph from the DataFrame
    G = nx.from_pandas_edgelist(bonds_df, 0, 1)  # Assuming columns 0 and 1 are the node pairings
    
    # Perform a BFS starting from the to_node, but avoiding from_node
    connected = set()
    queue = [to_node]
    visited = set([from_node])  # Start with the from_node to avoid it in the traversal
    
    while queue:
        current_node = queue.pop(0)
        if current_node in visited:
            continue
        visited.add(current_node)
        
        for neighbor in G.neighbors(current_node):
            if neighbor not in visited:
                connected.add(neighbor)
                queue.append(neighbor)
     
    # Return the set of connected nodes excluding from_node itself
    print("connected")
    print(connected)
    return connected
   
def list_xyz_files(directory):
    #List all .xyz files in the given directory.
    return [f for f in os.listdir(directory) if f.endswith('.xyz')]

def choose_file_menu(files):
    #Create a simple Tkinter dialog to choose a file from the list.
    root = Tk()
    root.withdraw()  # Hide the root window
    root.attributes('-topmost', True)  # Bring the dialog to the front

    # Create the menu string
    menu_str = "\n".join([f"{i}: {file}" for i, file in enumerate(files)])
    
# Function to open file dialog and select a file
def select_file():
    directory = r"/home/nati/Documents/GitHub/Alpha-Chemistry-UI/Roni Alpha/extrctfolder/Optimized_structures_xyz"

    root = Tk()
    root.withdraw()  # Hide the root window
    root.attributes('-topmost', True)  # Bring the file dialog to the front

    file_path = filedialog.askopenfilename(
        initialdir=directory,
        title="Select a molecule file",
        filetypes=(("XYZ files", ".xyz"), ("All files", ".*"))
    )
    root.destroy()  # Close the Tkinter root window

    # Get the index of the selected file
    if file_path:
        files = [f for f in os.listdir(directory) if f.endswith('.xyz') and os.path.isfile(os.path.join(directory, f))]
        print(files)
        file_name = os.path.basename(file_path)
        file_index = files.index(file_name) if file_name in files else -1
        return file_path, file_index
    else:
        return None, -1

     #Read data from the CSV file and store it in a dictionary
def read_element_data_from_csv(csv_file_path):
    element_data = {}
    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # Skip header row
        for row in csv_reader:
            atomic_number = int(row[0])
            color = row[2]
            radius = float(row[3])
            element_data[atomic_number] = {'color': color, 'radius': radius}
    return element_data
# my calculation of B1
def my_calc_B1(transformed_plane, avs, edited_coordinates_df, column_index):
    """
    Parameters
    ----------
    transformed_plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
            [-0.7384 -0.5135]
            [-0.3759 -0.271 ]
            [-1.1046 -0.8966]
            [ 0.6763  0.5885]
    avs : list
        the max & min of the [x,z] columns from the transformed_plane.
        example:[0.6763, -1.1046, 0.5885, -0.8966
                 ]
    edited_coordinates_df : TYPE
        DESCRIPTION.
    column_index : int
        0 or 1 depending- being used for transformed plane.
    """
    ## get the index of the min value in the column compared to the avs.min
    idx = np.where(np.isclose(np.abs(transformed_plane[:, column_index]), (avs.min()).round(4)))[0][0]
    if transformed_plane[idx, column_index] < 0:
        idx = np.where(np.isclose(transformed_plane[:, column_index], transformed_plane[:, column_index].min()))[0][0]
        bool_list = np.logical_and(transformed_plane[:, column_index] >= transformed_plane[idx, column_index],
                                   transformed_plane[:, column_index] <= transformed_plane[idx, column_index] + 1)

        transformed_plane[:, column_index] = -transformed_plane[:, column_index]
    else:
        bool_list = np.logical_and(transformed_plane[:, column_index] >= transformed_plane[idx, column_index] - 1,
                                   transformed_plane[:, column_index] <= transformed_plane[idx, column_index])

    against, against_loc = [], []
    B1, B1_loc = [], []
    B1_xz = [0,0 ]
    for i in range(1, transformed_plane.shape[0]):
        if bool_list[i]:
            against.append(np.array(transformed_plane[i, column_index] + edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        if len(against) > 0:
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
        else:
            B1.append(np.abs(transformed_plane[idx, column_index] + edited_coordinates_df['radius'].iloc[idx]))
            B1_loc.append(edited_coordinates_df['radius'].iloc[idx])

    
    
    return [B1, B1_loc, transformed_plane[idx]]


def my_b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane):
    """
    extended_df : pd.DataFrame
    b1s : list
    b1s_loc : list
    degree_list : list
    plane : np.array
    """
    degree = []
    b1s_xz = []
    for degree in degree_list:
        transformed_plane = get_transfomed_plane_for_sterimol(plane, degree)

        avs = np.abs([max(transformed_plane[:, 0]), min(transformed_plane[:, 0]),
                      max(transformed_plane[:, 1]), min(transformed_plane[:, 1])])

       # print("Transformed Plane:", transformed_plane)
        if min(avs) == 0:
            min_avs_indices = np.where(avs == min(avs))[0]
            if any(index in [0, 1] for index in min_avs_indices):
                tc = np.round(transformed_plane, 1)
                B1 = max(extended_df['radius'].iloc[np.where(tc[:, 0] == 0)])
                B1_loc = extended_df['L'].iloc[np.argmax(extended_df['radius'].iloc[np.where(tc[:, 0] == 0)])]
                b1s.append(B1)
                b1s_loc.append(B1_loc)
                b1s_xz.append(avs)
                continue  # Skip the rest of the loop

            elif any(index in [2, 3] for index in min_avs_indices):
                tc = np.round(transformed_plane, 1)
                B1 = max(extended_df['radius'].iloc[np.where(tc[:, 1] == 0)])
                B1_loc = extended_df['L'].iloc[np.argmax(extended_df['radius'].iloc[np.where(tc[:, 1] == 0)])]
                b1s.append(B1)
                b1s_loc.append(B1_loc)
                b1s_xz.append(avs)
                continue

        if np.where(avs == avs.min())[0][0] in [0, 1]:
            B1, B1_loc, B1_xz = my_calc_B1(transformed_plane, avs, extended_df, 0)

        elif np.where(avs == avs.min())[0][0] in [2, 3]:
            B1, B1_loc, B1_xz = my_calc_B1(transformed_plane, avs, extended_df, 1)

        b1s.append(np.unique(np.vstack(B1)).max())  ####check
        b1s_loc.append(np.unique(np.vstack(B1_loc)).max())
        b1s_xz.append(B1_xz)

    return b1s, b1s_loc, b1s_xz

#a function to get b1s
def my_get_b1s_list(extended_df, scans=90 // 5):
    b1s, b1s_loc = [], []
    scans = scans
    degree_list = list(range(18, 108, scans))
    plane = np.array(extended_df[['x', 'z']].astype(float))
    b1, b1_loc, b1_xz = my_b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
    print(f"Original plane vectors before transformation: {plane}")
    
    if b1s:
        try:
            back_ang = degree_list[np.where(b1s == min(b1s))[0][0]] - scans
            front_ang = degree_list[np.where(b1s == min(b1s))[0][0]] + scans
            degree_list = range(back_ang, front_ang + 1)
        except:
            print(np.where(np.isclose(b1s, min(b1s), atol=1e-8)))
            back_ang = degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]] - scans
            front_ang = degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]] + scans
            degree_list = range(back_ang, front_ang + 1)
    else:
        print('no b1s found')
        return np.array(b1s), np.array(b1s_loc), []

    # print(f'specific degree list: {degree_list}')
    b1, b1_loc, b1_xz = my_b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
    # print(f'b1 arrays: {[np.array(b1s),np.array(b1s_loc)]}')
    return np.array(b1s), np.array(b1s_loc), b1_xz
    # b1_vectors[b1s.index(min(b1s[b1s>=0]))]

# Create a dictionary to map atom names to atomic numbers
atom_names_to_numbers = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
    'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
    'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
    'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
    'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
    'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
    'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
    'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103
}

def filter_bonds(bonds_df, atom1_index, atom2_index):
    def remove_connections(bonds_df, atom1, atom2):
        # Find all connections of atom1
        connections = bonds_df[(bonds_df[0] == atom1) | (bonds_df[1] == atom1)]

        # Remove all connections of atom1 except for the connection with atom2
        bonds_df = bonds_df[~(((bonds_df[0] == atom1) | (bonds_df[1] == atom1)) &
                              (~((bonds_df[0] == atom2) | (bonds_df[1] == atom2))))]

        # Recursively remove connections of the connected atoms
        for _, row in connections.iterrows():
            connected_atom = row[1] if row[0] == atom1 else row[0]
            if connected_atom != atom2:
                bonds_df = remove_connections(bonds_df, connected_atom, atom1)

        return bonds_df

    return remove_connections(bonds_df, atom1_index, atom2_index)

def gather_indices(bonds_df):
    return list(set(map(int, bonds_df.values.flatten())))



# Function to plot molecule
def plot_molecule(xyz_data, connections, element_data, atom_numbers, file_path, file_name, nav):
    global real_xyz_data
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  # Create 3D axes

    atoms_plot_data = []
    
    cur_molecule = Molecule(file_path)
    
    # Function to handle mouse click event
    def onpick(event):
        global real_xyz_data
        if len(nav.selected_atoms) >= nav.num_atoms_to_pick:
            # Limit selection to the specified number of atoms
            return

        clicked_scatter = event.artist
        clicked_index = atoms_plot_data.index(clicked_scatter)
        chosen_atom_position = clicked_index + 1  # Index starts from 0, so add 1 for atom position
        nav.selected_atoms.append(chosen_atom_position)
        print(f"Selected atom position in file: {chosen_atom_position} ")
        
        # If two atoms are selected, draw the axis line and calculate Sterimol parameters
        if len(nav.selected_atoms) == 2:
            atom1_index = nav.selected_atoms[0] - 1  # Convert to zero-based index
            atom2_index = nav.selected_atoms[1] - 1  # Convert to zero-based index
            atom1_coords = real_xyz_data[atom1_index]
            atom2_coords = real_xyz_data[atom2_index]
            direction_vector = atom2_coords - atom1_coords
            direction_vector /= np.linalg.norm(direction_vector)  # Normalize the vector

            coords_df  = get_df_from_file(file_path)
            

            # Find the furthest atom in the direction of the second atom
            max_projection = max(np.dot((real_xyz_data - atom1_coords), direction_vector))
            extended_point = atom1_coords + max_projection * direction_vector            

            # Plot the line from the first atom to the furthest atom in the direction of the second atom
            ax.plot([atom1_coords[0], extended_point[0]],
                    [atom1_coords[1], extended_point[1]],
                    [atom1_coords[2], extended_point[2]], color='blue')

            # Add an arrowhead
            ax.quiver(atom2_coords[0], atom2_coords[1], atom2_coords[2],
                      direction_vector[0], direction_vector[1], direction_vector[2],
                      length=0.1, color='blue', arrow_length_ratio=3)

           # Calculate Sterimol parameters using Alpha11
            os.chdir(r'/home/nati/Roni/Roni Alpha/Optimized_structures_xyz')
            mols = Alpha11.Molecules(os.getcwd())
            sterimol_params = mols.get_sterimol_dict([atom1_index+1, atom2_index+1])
            file_name_without_extension = file_name[0:file_name.rfind('.')]
            sterimol_param = sterimol_params[file_name_without_extension]

            print(sterimol_param)
            print(atom1_index)
            print(atom2_index)            

            # Visualize B1 and B5
            L = np.linalg.norm(atom2_coords - atom1_coords)
            L_vector = atom2_coords - atom1_coords
            midpoint = (atom1_coords + atom2_coords) / 2
            B1 = sterimol_param['B1'].iloc[0]
            B5 = sterimol_param['B5'].iloc[0]

            # Filter atoms in the direction of the second atom
            bonds_df = cur_molecule.bonds_df
            """
            ask itay promisin to delete
            base_atoms = get_sterimol_indices(coords_df, bonds_df)
            bonds_direction = direction_atoms_for_sterimol(bonds_df, base_atoms)
            new_coordinates_df = preform_coordination_transformation(coords_df, bonds_direction)
            if True:
                connected_from_direction = get_molecule_connections(bonds_df, base_atoms[atom1_index], base_atoms[atom2_index])
            else:
                connected_from_direction = None
            bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, new_coordinates_df)
            
            extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, 'blue')
            edited_coords = filter_atoms_for_sterimol(bonded_atoms_df, coords_df)
            """
            
            indices = gather_indices(filter_bonds(bonds_df, 2, 1))
            indices_zero_based = [idx - 1 for idx in indices]
            edited_coordinates = coords_df.iloc[indices_zero_based]
            edited_coordinates['Original_Index'] = indices
            

            
            A = np.array([atom1_coords[0], atom1_coords[1], atom1_coords[2]])
            B = np.array([atom2_coords[0], atom2_coords[1], atom2_coords[2]])

            """
            ask itay
            # Function to project a point onto the line
            def get_project_magniute(P, A, B):
                            
                a = np.array(A - A)
                b = np.array(B - A)
                p = np.array(P - A)
                
                distance_b = np.linalg.norm(b)
                distance_p = np.linalg.norm(p)

                # dot between A and B
                above = p[0] * b[0] +  p[1] * b[1] + p[2] * b[2]
                lower = distance_b * distance_p
                alpha = np.arccos(above/lower)
                
                final_dis = distance_p * np.sin(alpha)
                
                return final_dis
            
            # Project each point and calculate the magnitudes of the projections
            coords_df['Projection Magnitude'] = coords_df.apply(
                lambda row: get_project_magniute(np.array([row['x'], row['y'], row['z']]), A, B), axis=1
            )
            """
            
            # Function to project a point onto the line
            def get_project_magniute(P, A, B):

                AP = P - A
                AB = B - A
                alpha = np.arccos(np.dot(AP, AB) / (np.linalg.norm(AP) * np.linalg.norm(AB)) )
                
                dis = np.linalg.norm(AP) * np.sin(alpha)                
                
                return dis

            def get_projection_vector(P, A, B):
                
                AP = P - A
                AB = B - A
                alpha = np.arccos(np.dot(AP, AB) / (np.linalg.norm(AP) * np.linalg.norm(AB)) )
                
                dis = np.linalg.norm(AP) * np.sin(alpha)                
                
                H = A + dis * AB
                
                return H - P


                    # Calculate a perpendicular vector to the blue axis
            def get_perpendicular_vector(direction_vector):
                if not np.allclose(direction_vector, [0, 0, 1]):
                    perp_vector = np.cross(direction_vector, [0, 0, 1])
                else:
                    perp_vector = np.cross(direction_vector, [0, 1, 0])
                perp_vector /= np.linalg.norm(perp_vector)
                return perp_vector

            perpendicular_vector = get_perpendicular_vector(direction_vector)
            perpendicular_vector *= B1 / 2  # Scale by B1/2 for visualization
           
    
            # calculate things to prepare for b1 Calculation
            base_atoms = get_sterimol_indices(coords_df, bonds_df)

            bonds_direction = direction_atoms_for_sterimol(bonds_df, base_atoms)
            new_coordinates_df = preform_coordination_transformation(coords_df, bonds_direction)
            connected_from_direction = get_molecule_connections(bonds_df, base_atoms[0], base_atoms[1])
            bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction,
                                                           new_coordinates_df)  ## xyz data frame with only the connected atoms

            extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, 'blue')
            edited_coords = filter_atoms_for_sterimol(bonded_atoms_df, coords_df)

            # calculate b1 vector
            edited_coordinates_df = filter_atoms_for_sterimol(bonded_atoms_df, extended_df)
            print("Edited", edited_coordinates_df)
            base_atoms = get_sterimol_indices(coords_df, bonds_df)
            print(sterimol_param['B1'])
            b1s,b1s_loc,b1s_xz=my_get_b1s_list(edited_coordinates_df)
            b1_index=np.where(b1s == min(b1s[b1s>=0]))[0][0]
            b1_xz = b1s_xz[b1_index]
            b1 = b1s[b1_index]
            b1_loc = b1s_loc[b1_index]
            b1_vector = np.array([b1_xz[0], b1_loc, b1_xz[1]])
            print("b1:", b1, "b1_xz:", b1_xz, "b1_loc:", b1_loc)
            
            coords_df = get_df_from_file(file_path)
            bonds_df = cur_molecule.bonds_df
            print("bonds_df")
            print(type(bonds_df))
            print(bonds_df)
            filter_indices=list(get_connected_nodes(bonds_df ,atom1_index+1 ,atom2_index+1))+[atom1_index+1 ,atom2_index+1]
            filter_indices = [num - 1 for num in filter_indices]
            filtered_df = edited_coordinates_df.loc[filter_indices]
            
            
            for _, row in filtered_df.iterrows():
                # Extract x, y, z values
                x, y, z = row['x'], row['y'], row['z']
                # Scatter each point with a yellow color and a radius of 2
                ax.scatter(x, y, z, color='yellow', s=100)
                
            
            
            edited_coordinates_df['Projection Magnitude'] = edited_coordinates_df.apply(
                lambda row: get_project_magniute(np.array([row['x'], row['y'], row['z']]), A, B), axis=1
            )

            max_projection_point = edited_coordinates_df.loc[edited_coordinates_df['Projection Magnitude'].idxmax()]

            ax.scatter(max_projection_point['x'], max_projection_point['y'], max_projection_point['z'], c='green', s=(radius * 100), picker=True)
            
            starting_point = (atom1_coords + atom2_coords) / 2  # Midpoint on the blue axis line
            ending_point = np.array([max_projection_point['x'], max_projection_point['y'],
                                     max_projection_point['z']])  # Coordinates of the atom with the radius of 100
            b5_vector = ending_point - starting_point
            L_vector_normalized = L_vector / np.linalg.norm(L_vector)
            projection = np.dot(b5_vector, L_vector_normalized) * L_vector_normalized
            new_b5_vector = b5_vector - projection
            new_ending_point = starting_point + new_b5_vector
            B5_loc = sterimol_param['loc_B5'].iloc[0]
            start_point_B5 = atom1_coords + B5_loc * direction_vector
            
            # Filter points where Projection Magnitude is greater than 0
            points_above_zero = edited_coordinates_df[edited_coordinates_df['Projection Magnitude'] > 0]

            # Plot all points where Projection Magnitude is greater than 0
            ax.scatter(
                points_above_zero['x'], 
                points_above_zero['y'], 
                points_above_zero['z'], 
                c='blue',  # Choose a color for these points
                s=20,      # Set the size of each point
                picker=True
            )

            # Define a perpendicular direction vector for B1 line
            # Assuming 'perpendicular_direction' is orthogonal to 'direction_vector'
            # (you may already have this calculation in place; otherwise, compute it)
            perpendicular_direction = np.array([-direction_vector[1], direction_vector[0], 0])
            perpendicular_direction /= np.linalg.norm(perpendicular_direction)

            # Calculate the end point of the B1 line
            end_point_B5 = start_point_B5 + b5_vector

            # Draw the B1 line in red
            ax.plot([start_point_B5[0], end_point_B5[0]],
                    [start_point_B5[1], end_point_B5[1]],
                    [start_point_B5[2], end_point_B5[2]], color='green', linestyle='--')
            
            
            
            
                # Define the blue axis (new y-axis) as the vector between atom1_coords and atom2_coords
            L_vector = atom2_coords - atom1_coords  # Blue axis
            L_vector_normalized = L_vector / np.linalg.norm(L_vector)  # Normalize to make it a unit vector

            # Project b1_loc onto the blue axis to find the starting point on the y-axis
            b1_loc_vector = np.array([b1_xz[0], b1_loc, b1_xz[1]])  # Position of the B1 vector

            # Find the projection of b1_loc onto the blue axis (y-axis)
            projection_on_y_axis = atom1_coords + np.dot(b1_loc_vector - atom1_coords, L_vector_normalized) * L_vector_normalized

            # Now calculate the perpendicular component of the B1 vector relative to the blue axis
            projection = np.dot(b1_vector, L_vector_normalized) * L_vector_normalized  # Part of B1 aligned with blue axis
            perpendicular_b1_vector = b1_vector - projection  # Perpendicular part

            # Plot the B1 vector, starting from the projection on the y-axis (along the blue line)
            start_point = projection_on_y_axis  # This is now aligned with the blue axis
            end_point = start_point + perpendicular_b1_vector  # End point of the B1 vector
           
            # Project b1_vector perpendicular to the blue line
            projection = np.dot(b1_vector, L_vector_normalized) * L_vector_normalized
            perpendicular_b1_vector = b1_vector - projection
            b1_horrible_vector = np.array([b1_xz[0], b1_loc, b1_xz[1]])
            b1_normalized_vector = b1_horrible_vector / np.sqrt(np.sum(b1_horrible_vector**2))
            b1_vector = b1_normalized_vector * b1
            B1_loc = sterimol_param['loc_B1'].iloc[0]
           
            # Calculate the starting point along the blue line at B1_loc distance
            start_point_B1 = atom1_coords + B1_loc * direction_vector

            # Define a perpendicular direction vector for B1 line
            # Assuming 'perpendicular_direction' is orthogonal to 'direction_vector'
            # (you may already have this calculation in place; otherwise, compute it)
            perpendicular_direction = np.array([-direction_vector[1], direction_vector[0], 0])
            perpendicular_direction /= np.linalg.norm(perpendicular_direction)

            # Calculate the end point of the B1 line
            end_point_B1 = start_point_B1 + b1_vector

            # Draw the B1 line in red
            ax.plot([start_point_B1[0], end_point_B1[0]],
                    [start_point_B1[1], end_point_B1[1]],
                    [start_point_B1[2], end_point_B1[2]], color='red', linestyle='--')
            
                    
            # Optionally, add an arrowhead for the purple line
            ax.quiver(start_point[0], start_point[1], start_point[2],
                      perpendicular_vector[0], perpendicular_vector[1], perpendicular_vector[2],
                      length=0.1, color='red', arrow_length_ratio=0.2)

           

        # If all atoms are selected, display the list of selected atom positions
        if len(nav.selected_atoms) == nav.num_atoms_to_pick:
            selected_positions = ", ".join(map(str, nav.selected_atoms))
            nav.save_selected_atoms(selected_positions)  # Save selected atom positions to a file
            messagebox.showinfo("Selected Atom Positions", f"You have selected the following atom positions: {selected_positions}")
            # Clear the selected atoms list for the next molecule
            nav.selected_atoms.clear()

    # Scatter plot for atoms
    coords_df1  = get_df_from_file(file_path)  
    
    bonds_df1 = cur_molecule.bonds_df

    base_atoms1 = get_sterimol_indices(coords_df1, bonds_df1)
    bonds_direction1 = direction_atoms_for_sterimol(bonds_df1, base_atoms1)
    new_coordinates_df1 = preform_coordination_transformation(coords_df1, bonds_direction1)

    connected_from_direction1 = get_molecule_connections(bonds_df1, base_atoms1[0], base_atoms1[1])
    bonded_atoms_df1 = get_specific_bonded_atoms_df(bonds_df1, connected_from_direction1, new_coordinates_df1)

    extended_df1 = get_extended_df_for_sterimol(new_coordinates_df1, bonds_df1, 'blue')
    edited_coords1 = filter_atoms_for_sterimol(bonded_atoms_df1, coords_df1)

    # Edited coordinates
    edited_coordinates_df1 = filter_atoms_for_sterimol(bonded_atoms_df1, extended_df1)
    
    # Scatter plot for atoms
    print("world xyz")
    print(xyz_data)
    print(type(xyz_data))
    xyz_data = edited_coordinates_df1[['x', 'y', 'z']].to_numpy()
    real_xyz_data = xyz_data

   # real_xyz_data[[['x','y','z']]] = new_coordinates_df

    print(real_xyz_data)
    for i, coords in enumerate(xyz_data):
        element_info = element_data.get(atom_numbers[i])
        color = element_info['color'] if element_info else 'blue'
        radius = element_info['radius'] if element_info else 1.0
        atoms_plot_data.append(ax.scatter(coords[0], coords[1], coords[2], c=color, s=radius ** 4, picker=True))

    # Plot connections between atoms
    for connection in connections:
        atom1_index, atom2_index = connection
        atom1_coords = xyz_data[atom1_index]
        atom2_coords = xyz_data[atom2_index]
        ax.plot([atom1_coords[0], atom2_coords[0]],
                [atom1_coords[1], atom2_coords[1]],
                [atom1_coords[2], atom2_coords[2]], color='black')

    ax.grid(False)  # Disable grid lines
    ax.set_aspect('equal')  # Set aspect ratio to equal

    # Connect click event to callback function
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.title(file_name)
    ax.axis('off')  # Turn off axes
    plt.show()

# Rest of the code remains the same
# Function to remove violating connections
def remove_violating_connections(xyz_data, connections, atom_symbols, threshold_distance):
    # Create a list to store the filtered connections
    filtered_connections = []
    for connection in connections:
        atom1_index, atom2_index = connection
        atom1_symbol = atom_symbols[atom1_index]
        atom2_symbol = atom_symbols[atom2_index]
        # Check if the connection should be kept based on the atom symbols
        if not ((atom1_symbol == 'H' and atom2_symbol not in ['N', 'O', 'F']) or
                (atom2_symbol == 'H' and atom1_symbol not in ['N', 'O', 'F']) or
                (atom1_symbol == 'H' and atom2_symbol == 'H')):
            filtered_connections.append(connection)
    return filtered_connections

# Class to manage navigation
# Class to manage navigation
class Navigation:
    def __init__(self, all_files, extracted_folder, element_data, num_atoms_to_pick, index=0):
        self.all_files = all_files
        self.extracted_folder = extracted_folder
        self.element_data = element_data
        self.current_index = index
        self.num_files = len(all_files)
        #self.fig, self.ax = plt.subplots()
       # self.axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
       # self.axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
       # self.bnext = Button(self.axnext, 'Next')
       # self.bprev = Button(self.axprev, 'Previous')
       # self.bnext.on_clicked(self.next_molecule)
       # self.bprev.on_clicked(self.prev_molecule)
        self.num_atoms_to_pick = num_atoms_to_pick
        self.selected_atoms = []  # Initialize the list of selected atoms
        self.list_number = 1  # Initialize list number
        self.plot_current_molecule()

    def next_molecule(self, event):
        self.current_index = (self.current_index + 1) % self.num_files
        self.plot_current_molecule()

    def prev_molecule(self, event):
        self.current_index = (self.current_index - 1) % self.num_files
        self.plot_current_molecule()

    def plot_current_molecule(self):
        file_name = self.all_files[self.current_index]
        file_path = os.path.join(self.extracted_folder, file_name)

        # Read XYZ data
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Assuming the format is "Element X Y Z"
        xyz_data = []
        atom_numbers = []  # Store atom numbers
        atom_symbols = []  # Store atom symbols
        for i, line in enumerate(lines, 1):
            data = line.strip().split()
            if len(data) >= 4:  # Ensure each line has at least four elements
                atom_name = data[0]  # Atom name
                atomic_number = atom_names_to_numbers.get(atom_name, 6)  # Default atomic number to 6 if not found
                coordinates = list(map(float, data[1:]))  # Convert remaining elements to float
                xyz_data.append(coordinates)
                atom_numbers.append(atomic_number)  # Store atomic number
                atom_symbols.append(atom_name)  # Store atom symbol

        xyz_data = np.array(xyz_data)

        # Compute distances between points
        connections = []
        for pair in combinations(range(len(xyz_data)), 2):
            dist = np.linalg.norm(xyz_data[pair[0]] - xyz_data[pair[1]])
            if dist < 1.8:
                connections.append((pair[0], pair[1]))

        # Remove violating connections
        connections = remove_violating_connections(xyz_data, connections, atom_symbols, threshold_distance=2.12)

        plot_molecule(xyz_data, connections, self.element_data, atom_numbers, file_path, file_name, self)  # Pass self as nav

    def save_selected_atoms(self, selected_positions):
        with open('selected_atom_positions.txt', 'a') as f:  # Use 'a' mode to append
            f.write(f"List {self.list_number}:\n")  # Add a header for the new list
            f.write(selected_positions + '\n')  # Append selected positions
        self.list_number += 1  # Increment list number for the next list

# Function to ask for the number of atoms to pick using a GUI window
def ask_num_atoms():
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    num_atoms = simpledialog.askinteger("Number of Atoms", "Enter the number of atoms you want to pick:")
    return num_atoms



# Main function
def main ():
    # Path to the CSV file containing element data
    csv_file_path ='/home/nati/Documents/GitHub/Alpha-Chemistry-UI/Roni Alpha/rsdii_table.csv'
    # Read element data from CSV
    element_data = read_element_data_from_csv(csv_file_path)

    # Path to the zip file containing XYZ data
    zip_file_path = r"/home/nati/Documents/GitHub/Alpha-Chemistry-UI/Roni Alpha/Optimized_structures_xyz.zip"
    extracted_folder = r"/home/nati/Documents/GitHub/Alpha-Chemistry-UI/Roni Alpha/extrctfolder/Optimized_structures_xyz"

    # Extract all files from the zip
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(extracted_folder)

    all_files = [o for o in os.listdir(extracted_folder) if o[-4:]==".xyz"]

    # Select the file
    selected_file, selected_index = select_file()
    print(selected_file)
    directory = r"directory = r'/home/nati/Documents/GitHub/Alpha-Chemistry-UI/Roni Alpha/extrctfolder/Optimized_structures_xyz"
    print(selected_index)
    # Ask for the number of atoms to pick
    num_atoms_to_pick = ask_num_atoms()

    nav = Navigation(all_files, extracted_folder, element_data, num_atoms_to_pick, index=selected_index)


    while True:
        if messagebox.askyesno("Continue", "Do you want to pick atoms in another molecule?"):
            num_atoms_to_pick = ask_num_atoms()  # Ask for the number of atoms to pick
            nav.num_atoms_to_pick = num_atoms_to_pick  # Update the number of atoms to pick
            nav.selected_atoms = []  # Reset the list of selected atoms
            nav.plot_current_molecule()  # Start the selection process
        else:
            break
            
# Entry point of the script
if __name__ == "__main__":
    main()