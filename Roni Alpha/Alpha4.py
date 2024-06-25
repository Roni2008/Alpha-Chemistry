import zipfile
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os
import csv
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Button
from mpl_toolkits.mplot3d import proj3d  # Import proj3d explicitly

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

# Function to plot molecule
def plot_molecule(xyz_data, connections, element_data, file_name, atom_numbers):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  # Create 3D axes

    atoms_plot_data = []

    # Function to handle mouse click event
    def onpick(event):
        clicked_scatter = event.artist
        for item in atoms_plot_data:
            if item == clicked_scatter:
                clicked_index = atoms_plot_data.index(item)
        chosen_atom = atom_numbers[clicked_index]

        # Open a window displaying the atom number
        plt.figure()
        plt.text(0.5, 0.5, str(chosen_atom), fontsize=20, ha='center')
        plt.axis('off')
        plt.show()

    # Scatter plot for atoms
    for i, coords in enumerate(xyz_data):
        element_info = element_data.get(atom_numbers[i])
        color = element_info['color'] if element_info else 'blue'
        radius = element_info['radius'] if element_info else 1.0
        atoms_plot_data.append(ax.scatter(coords[0], coords[1], coords[2], c=color, s=radius ** 2,picker=True, label=str(atom_numbers[i])))

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(file_name)  # Set title to file name

    # Plot connections between atoms
    for connection in connections:
        atom1_index, atom2_index = connection
        atom1_coords = xyz_data[atom1_index - 1]  # Adjust index to start from 0
        atom2_coords = xyz_data[atom2_index - 1]  # Adjust index to start from 0
        ax.plot([atom1_coords[0], atom2_coords[0]],
                [atom1_coords[1], atom2_coords[1]],
                [atom1_coords[2], atom2_coords[2]], color='black')

    ax.frame_on = False  # Hide axes lines and box
    fig.canvas.mpl_connect('pick_event', onpick)  # Connect click event to callback function
    plt.show()

# Class to manage navigation
class Navigation:
    def __init__(self, all_files, extracted_folder, element_data):
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
        atom_numbers = []
        for i, line in enumerate(lines, 1):
            data = line.strip().split()
            if len(data) >= 4:  # Ensure each line has at least four elements
                element = data[0]  
                coordinates = list(map(float, data[1:]))  # Convert remaining elements to float
                xyz_data.append(coordinates)
                atom_numbers.append(i)  # Store atom number

        xyz_data = np.array(xyz_data)

        # Compute distances between points
        connections = []
        for pair in combinations(range(len(xyz_data)), 2):
            if np.linalg.norm(xyz_data[pair[0]] - xyz_data[pair[1]]) <  1.8:
                connections.append(pair)

        plot_molecule(xyz_data, connections, self.element_data, file_name, atom_numbers)

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
    nav = Navigation(all_files, extracted_folder, element_data)

# Entry point of the script
if __name__ == "__main__":
    main()