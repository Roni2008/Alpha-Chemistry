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

# Itay Added A comment!!!!! (Example)

# Read data from the CSV file and store it in a dictionary
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

# Function to plot molecule
def plot_molecule(xyz_data, connections, element_data, atom_numbers, file_name, nav):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  # Create 3D axes

    atoms_plot_data = []

    # Function to handle mouse click event
    def onpick(event):
        if len(nav.selected_atoms) >= nav.num_atoms_to_pick:
            # Limit selection to the specified number of atoms
            return
    
        clicked_scatter = event.artist
        clicked_index = atoms_plot_data.index(clicked_scatter)
        chosen_atom_position = clicked_index + 1  # Index starts from 0, so add 1 for atom position
        nav.selected_atoms.append(chosen_atom_position)
        print(f"Selected atom position in file: {chosen_atom_position}")
        
         

        # If two atoms are selected, draw the axis line and calculate Sterimol parameters
        if len(nav.selected_atoms) == 2:
            atom1_index = nav.selected_atoms[0] - 1  # Convert to zero-based index
            atom2_index = nav.selected_atoms[1] - 1  # Convert to zero-based index
            atom1_coords = xyz_data[atom1_index]
            atom2_coords = xyz_data[atom2_index]
            direction_vector = atom2_coords - atom1_coords
            direction_vector /= np.linalg.norm(direction_vector)  # Normalize the vector
            
             # Find the furthest atom in the direction of the second atom
            max_projection = max(np.dot((xyz_data - atom1_coords), direction_vector))
            extended_point = atom1_coords + max_projection * direction_vector
            
             # Plot the line from the first atom to the furthest atom in the direction of the second atom
            ax.plot([atom1_coords[0], extended_point[0]],
                    [atom1_coords[1], extended_point[1]],
                    [atom1_coords[2], extended_point[2]], color='blue')

            # Add an arrowhead
            ax.quiver(atom2_coords[0], atom2_coords[1], atom2_coords[2], 
                      direction_vector[0], direction_vector[1], direction_vector[2], 
                      length=0.1, color='blue', arrow_length_ratio=3)
        

            """
            L = np.linalg.norm(atom2_coords - atom1_coords)
            distances = np.linalg.norm(xyz_data - (atom1_coords + atom2_coords) / 2, axis=1)
            B1 = min(distances)
            B5 = max(distances)
            sterimol_params = f"L: {L:.2f}, B1: {B1:.2f}, B5: {B5:.2f}"
            print(f"Sterimol Parameters: {sterimol_params}")
            messagebox.showinfo("Sterimol Parameters", f"Sterimol Parameters: {sterimol_params}")
            plt.draw()
            """

            # Calculate Sterimol parameters
            os.chdir(r'/home/nati/Roni/Roni Alpha/Optimized_structures_xyz')
            mols=Alpha11.Molecules(os.getcwd())
            sterimol_params = mols.get_sterimol_dict([4,7])
            file_name_without_extension = file_name[0:file_name.rfind('.')]
            sterimol_param = sterimol_params[file_name_without_extension]
            print(sterimol_param)
            
            
           
            # Visualize B1 and B5
            L = np.linalg.norm(atom2_coords - atom1_coords)
            midpoint = (atom1_coords + atom2_coords) / 2
            B1 = sterimol_param['B1'].iloc[0]
            B5 = sterimol_param['B5'].iloc[0]

            # Generate two perpendicular vectors in the plane perpendicular to the direction_vector
            random_vector = np.random.rand(3)  # Generate a random vector
            perp_vector1 = np.cross(direction_vector, random_vector)
            perp_vector1 /= np.linalg.norm(perp_vector1)  # Normalize the vector
            perp_vector2 = np.cross(direction_vector, perp_vector1)
            perp_vector2 /= np.linalg.norm(perp_vector2)  # Normalize the vector

            # Plot B1 vector in green
            ax.quiver(midpoint[0], midpoint[1], midpoint[2], 
                      perp_vector1[0], perp_vector1[1], perp_vector1[2], 
                      length=B1, color='green', arrow_length_ratio=0.1)

            # Plot B5 vector in red
            ax.quiver(midpoint[0], midpoint[1], midpoint[2], 
                      perp_vector2[0], perp_vector2[1], perp_vector2[2], 
                      length=B5, color='red', arrow_length_ratio=0.1)

            messagebox.showinfo("Sterimol Parameters", f"Sterimol Parameters: L: {L:.2f}, B1: {B1:.2f}, B5: {B5:.2f}")
            plt.draw()
            
            
            
        # If all atoms are selected, display the list of selected atom positions    
        if len(nav.selected_atoms) == nav.num_atoms_to_pick:
            selected_positions = ", ".join(map(str, nav.selected_atoms))
            nav.save_selected_atoms(selected_positions)  # Save selected atom positions to a file
            messagebox.showinfo("Selected Atom Positions", f"You have selected the following atom positions: {selected_positions}")
            # Clear the selected atoms list for the next molecule
            nav.selected_atoms.clear()

    # Scatter plot for atoms
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
class Navigation:
    def __init__(self, all_files, extracted_folder, element_data, num_atoms_to_pick):
        self.all_files = all_files
        self.extracted_folder = extracted_folder
        self.element_data = element_data
        self.current_index = 0
        self.num_files = len(all_files)
        self.fig, self.ax = plt.subplots()
        self.axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
        self.axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
        self.bnext = Button(self.axnext, 'Next')
        self.bprev = Button(self.axprev, 'Previous')
        self.bnext.on_clicked(self.next_molecule)
        self.bprev.on_clicked(self.prev_molecule)
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

        plot_molecule(xyz_data, connections, self.element_data, atom_numbers, file_name, self)  # Pass self as nav

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
def main():
    # Path to the CSV file containing element data
    csv_file_path = '/home/nati/Roni/Roni Alpha/rsdii_table.csv'
    # Read element data from CSV
    element_data = read_element_data_from_csv(csv_file_path)

    # Path to the zip file containing XYZ data
    zip_file_path = "/home/nati/Roni/Roni Alpha/Optimized_structures_xyz.zip"
    extracted_folder = "/home/nati/Roni/Roni Alpha/Optimized_structures_xyz"

    # Extract all files from the zip
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(extracted_folder)

    all_files = os.listdir(extracted_folder)

    # Ask for the number of atoms to pick
    num_atoms_to_pick = ask_num_atoms()

    nav = Navigation(all_files, extracted_folder, element_data, num_atoms_to_pick)

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
