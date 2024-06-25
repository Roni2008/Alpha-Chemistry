import zipfile
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os

# Define ELEMENT_DATA dictionary here
ELEMENT_DATA = {
    1: {'name': 'Hydrogen', 'color': '#ffffff', 'radius': 1.2},
    2: {'name': 'Helium', 'color': '#ffffff', 'radius': 1.43},
    3: {'name': 'Lithium', 'color': '#ffc0cb', 'radius': 2.12},
    4: {'name': 'Beryllium', 'color': '#b22222', 'radius': 1.98},
    5: {'name': 'Boron', 'color': '#ff1493', 'radius': 1.91},
    6: {'name': 'Carbon', 'color': '#00ff00', 'radius': 1.77},
    7: {'name': 'Nitrogen', 'color': '#c8c8c8', 'radius': 1.66},
    8: {'name': 'Oxygen', 'color': '#8f8fff', 'radius': 1.5},
    9: {'name': 'Fluorine', 'color': '#f00000', 'radius': 1.46},
    10: {'name': 'Neon', 'color': '#daa520', 'radius': 1.58},
    11: {'name': 'Sodium', 'color': '#ff1493', 'radius': 2.5},
    12: {'name': 'Magnesium', 'color': '#0000ff', 'radius': 2.51},
    13: {'name': 'Aluminum', 'color': '#228b22', 'radius': 2.25},
    14: {'name': 'Silicon', 'color': '#808090', 'radius': 2.19},
    15: {'name': 'Phosphorus', 'color': '#daa520', 'radius': 1.9},
    16: {'name': 'Sulfur', 'color': '#ffa500', 'radius': 1.89},
    17: {'name': 'Chlorine', 'color': '#ffc832', 'radius': 1.82},
    18: {'name': 'Argon', 'color': '#00ff00', 'radius': 1.94},
    19: {'name': 'Potassium', 'color': '#ff1493', 'radius': 2.73},
    20: {'name': 'Calcium', 'color': '#ff1493', 'radius': 2.62},
    # Continue with the rest of the elements...
    83: {'name': 'Bismuth', 'color': None, 'radius': None},
    84: {'name': 'Polonium', 'color': None, 'radius': None},
    85: {'name': 'Astatine', 'color': None, 'radius': None},
    87: {'name': 'Francium', 'color': None, 'radius': None},
    100: {'name': 'Fermium', 'color': None, 'radius': None},
    101: {'name': 'Mendelevium', 'color': None, 'radius': None},
    102: {'name': 'Nobelium', 'color': None, 'radius': None},
    103: {'name': 'Lawrencium', 'color': None, 'radius': None},
    104: {'name': 'Rutherfordium', 'color': None, 'radius': None},
    105: {'name': 'Dubnium', 'color': None, 'radius': None},
    106: {'name': 'Seaborgium', 'color': None, 'radius': None},
    107: {'name': 'Bohrium', 'color': None, 'radius': None},
    108: {'name': 'Hassium', 'color': None, 'radius': None},
    109: {'name': 'Meitnerium', 'color': None, 'radius': None},
    110: {'name': 'Darmstadtium', 'color': None, 'radius': None},
    111: {'name': 'Roentgenium', 'color': None, 'radius': None},
    112: {'name': 'Copernicium', 'color': None, 'radius': None},
    113: {'name': 'Nihonium', 'color': None, 'radius': None},
    114: {'name': 'Flerovium', 'color': None, 'radius': None},
    115: {'name': 'Moscovium', 'color': None, 'radius': None},
    116: {'name': 'Livermorium', 'color': None, 'radius': None},
    117: {'name': 'Tennessine', 'color': None, 'radius': None},
    118: {'name': 'Oganesson', 'color': None, 'radius': None}
}

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
        element_info = ELEMENT_DATA.get(i + 1)  # Element data starts from 1, not 0
        if element_info:
            color = element_info['color']
            radius = element_info['radius']
            ax.scatter(xyz_data[i, 0], xyz_data[i, 1], xyz_data[i, 2], c=color, s=radius**2, label=element_info['name'])
    for connection in connections:
        ax.plot([connection[0][0], connection[1][0]],
                [connection[0][1], connection[1][1]],
                [connection[0][2], connection[1][2]],
                color='black')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()

    input("Press Enter to view the next plot...")  
