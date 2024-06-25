import zipfile
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os

zip_file_path = "/home/nati/Roni/Roni Alpha/Optimized_structures_xyz.zip"
extracted_folder = "/home/nati/Roni/Roni Alpha/Optimized_structures_xyz"

# Extract all files from the zip
with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
    zip_ref.extractall(extracted_folder)

all_files = os.listdir(extracted_folder)
for file_name in all_files:
    file_path = os.path.join(extracted_folder, file_name)

    # Step 2: Read XYZ data
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Assuming the format is "Element X Y Z"
    elements = []
    xyz_data = []
    for line in lines:
        data = line.strip().split()
        if len(data) >= 4:  # Ensure each line has at least four elements
            element = data[0]  
            coordinates = list(map(float, data[1:]))  # Convert remaining elements to float
            elements.append(element)
            xyz_data.append(coordinates)

    xyz_data = np.array(xyz_data)

    # Step 3: Compute distances between points
    def distance(p1, p2):
        return np.linalg.norm(p2 - p1)

    connections = []
    for pair in combinations(xyz_data, 2):
        if distance(pair[0], pair[1]) < 1.8:
            connections.append(pair)

    connections = np.array(connections)

    # Step 4: Plot scatter graph and connect points
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(xyz_data)):
        ax.scatter(xyz_data[i, 0], xyz_data[i, 1], xyz_data[i, 2], c='b', marker='o')
        ax.text(xyz_data[i, 0], xyz_data[i, 1], xyz_data[i, 2], elements[i], color='red')
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