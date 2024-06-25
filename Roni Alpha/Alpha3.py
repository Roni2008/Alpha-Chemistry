import zipfile
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os
import csv

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

zip_file_path = "/home/nati/Roni/Roni Alpha/Optimized_structures_xyz.zip"
extracted_folder = "/home/nati/Roni/Roni Alpha/Optimized_structures_xyz"
csv_file_path = '/home/nati/Roni/Roni Alpha/rsdii_table.csv'

# Extract all files from the zip
with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
    zip_ref.extractall(extracted_folder)

# Read element data from CSV file
element_data = read_element_data_from_csv(csv_file_path)

all_files = os.listdir(extracted_folder)
for file_name in all_files:
    file_path = os.path.join(extracted_folder, file_name)

    # Read XYZ data
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Assuming the format is "Element X Y Z"
    xyz_data = []
    for line in lines:
        data = line.strip().split()
        if len(data) >= 4:  # Ensure each line has at least four elements
            coordinates = list(map(float, data[1:]))  # Convert remaining elements to float
            xyz_data.append(coordinates)

    xyz_data = np.array(xyz_data)

    # Compute distances between points
    def distance(p1, p2):
        return np.linalg.norm(p2 - p1)

    connections = []
    for pair in combinations(xyz_data, 2):
        if distance(pair[0], pair[1]) < 1.8:
            connections.append(pair)

    connections = np.array(connections)

    # Plot scatter graph and connect points
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(xyz_data)):
        color = 'blue'  # Default color
        radius = 1.0  # Default radius
        element_info = element_data.get(i + 1)  # Element data starts from 1, not 0
        if element_info:
            color = element_info['color']
            radius = element_info['radius']
        ax.scatter(xyz_data[i, 0], xyz_data[i, 1], xyz_data[i, 2], c=color, s=radius**2)
    for connection in connections:
        ax.plot([connection[0][0], connection[1][0]],
                [connection[0][1], connection[1][1]],
                [connection[0][2], connection[1][2]],
                color='black')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

    input("Press Enter to view the next plot...")

