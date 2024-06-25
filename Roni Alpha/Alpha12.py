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
import tempfile

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
                      length=0.1, color='blue', arrow_length_ratio=0.2)

            plt.draw()

            # Calculate and display Sterimol parameters
            calculate_sterimol_parameters(xyz_data, nav.selected_atoms, atom_numbers)

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
        atoms_plot_data.append(ax.scatter(coords[0], coords[1], coords[2], c=color, s=radius ** 2, picker=True))

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

# Calculate Sterimol parameters
def calculate_sterimol_parameters(xyz_data, selected_atoms, atom_numbers):
    atom1_index = selected_atoms[0] - 1
    atom2_index = selected_atoms[1] - 1
    atom1_coords = xyz_data[atom1_index]
    atom2_coords = xyz_data[atom2_index]
    direction_vector = atom2_coords - atom1_coords
    direction_vector /= np.linalg.norm(direction_vector)  # Normalize the vector

    # Calculate L
    L = np.linalg.norm(atom2_coords - atom1_coords)

    # Calculate B1 and B5
    projections = np.dot(xyz_data - atom1_coords, direction_vector)
    perpendicular_distances = np.linalg.norm(xyz_data - atom1_coords - np.outer(projections, direction_vector), axis=1)
    B1 = np.min(perpendicular_distances)
    B5 = np.max(perpendicular_distances)

    print(f"Sterimol Parameters:\nB1: {B1:.3f} Å\nB5: {B5:.3f} Å\nL: {L:.3f} Å")

# Function to remove violating connections
def remove_violating_connections(xyz_data, connections, atom_symbols, threshold_distance):
    filtered_connections = []
    for connection in connections:
        atom1_index, atom2_index = connection
        distance = np.linalg.norm(xyz_data[atom1_index] - xyz_data[atom2_index])
        if distance < threshold_distance:
            filtered_connections.append(connection)
    return filtered_connections

# Class to handle molecule navigation and atom selection
class MoleculeNavigator:
    def __init__(self, xyz_files, element_data):
        self.xyz_files = xyz_files
        self.element_data = element_data
        self.current_index = 0
        self.num_atoms_to_pick = 2  # Number of atoms to select
        self.selected_atoms = []

    def load_next_molecule(self):
        self.current_index = (self.current_index + 1) % len(self.xyz_files)
        self.display_current_molecule()

    def load_previous_molecule(self):
        self.current_index = (self.current_index - 1) % len(self.xyz_files)
        self.display_current_molecule()

    def display_current_molecule(self):
        xyz_file_path = self.xyz_files[self.current_index]
        xyz_data, atom_numbers, atom_symbols = read_xyz_file(xyz_file_path)
        connections = generate_connections(xyz_data, atom_symbols, threshold_distance=1.6)
        connections = remove_violating_connections(xyz_data, connections, atom_symbols, threshold_distance=4.0)
        file_name = os.path.basename(xyz_file_path)
        plot_molecule(xyz_data, connections, self.element_data, atom_numbers, file_name, self)

    def save_selected_atoms(self, selected_atoms):
        # Implement saving selected atom positions if needed
        print(f"Selected atom positions: {selected_atoms}")

def on_next(event, navigator):
    navigator.load_next_molecule()

def on_previous(event, navigator):
    navigator.load_previous_molecule()

# Function to read XYZ file
def read_xyz_file(xyz_file_path):
    with open(xyz_file_path, 'r') as file:  # Corrected unmatched parenthesis here
        lines = file.readlines()
    num_atoms = int(lines[0].strip())
    atom_data = lines[2:num_atoms + 2]
    xyz_data = []
    atom_numbers = []
    atom_symbols = []
    for atom_line in atom_data:
        parts = atom_line.split()
        atom_symbol = parts[0]
        x, y, z = map(float, parts[1:])
        xyz_data.append([x, y, z])
        atom_numbers.append(atom_names_to_numbers[atom_symbol])
        atom_symbols.append(atom_symbol)
    return np.array(xyz_data), atom_numbers, atom_symbols

# Function to generate connections
def generate_connections(xyz_data, atom_symbols, threshold_distance):
    connections = []
    for i, j in combinations(range(len(xyz_data)), 2):
        distance = np.linalg.norm(xyz_data[i] - xyz_data[j])
        if distance < threshold_distance:
            connections.append((i, j))
    return connections

# Main function to run the script
def main():
    csv_file_path = '/home/nati/Roni/Roni Alpha/rsdii_table.csv'
    zip_file_path = '/home/nati/Roni/Roni Alpha/Optimized_structures_xyz.zip'
    if not os.path.isfile(csv_file_path):
        print(f"CSV file not found: {csv_file_path}")
        return

    if not os.path.isfile(zip_file_path):
        print(f"ZIP file not found: {zip_file_path}")
        return

    element_data = read_element_data_from_csv(csv_file_path)

    # Extract the ZIP file to a temporary directory
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        with tempfile.TemporaryDirectory() as temp_dir:
            zip_ref.extractall(temp_dir)
            xyz_files = []
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    if file.endswith('.xyz'):
                        xyz_files.append(os.path.join(root, file))

            # Debugging statements
            print(f"Temporary directory: {temp_dir}")
            print(f"XYZ files found: {xyz_files}")

            if not xyz_files:
                print("No .xyz files found in the ZIP archive.")
                return  # Exit if no .xyz files are found

            navigator = MoleculeNavigator(xyz_files, element_data)

            # Create main window
            root = tk.Tk()
            root.title("Molecule Viewer")
            frame = tk.Frame(root)
            frame.pack()

            next_button = tk.Button(frame, text="Next", command=lambda: on_next(None, navigator))
            next_button.pack(side=tk.LEFT)

            previous_button = tk.Button(frame, text="Previous", command=lambda: on_previous(None, navigator))
            previous_button.pack(side=tk.LEFT)

            navigator.display_current_molecule()  # Display the first molecule

            root.mainloop()

if __name__ == "__main__":
    main()


