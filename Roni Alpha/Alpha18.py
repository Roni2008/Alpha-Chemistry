import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
import numpy as np

# Generate 30 random 3D points
np.random.seed(0)
points = np.random.rand(30, 3)

def transform_coordinates_b1_vector(points, idx1, idx2):
    """Transform points to align selected points with new axes."""
    # Translate so that idx1 is the origin
    origin = points[idx1]
    translated_points = points - origin

    # Vector from idx1 to idx2
    target_vector = translated_points[idx2]

    # Normalize the target vector
    target_vector = target_vector / np.linalg.norm(target_vector)

    # Construct the rotation matrix to align target_vector with Z-axis
    z_axis = np.array([0, 0, 1])
    v = np.cross(target_vector, z_axis)
    c = np.dot(target_vector, z_axis)
    s = np.linalg.norm(v)
    if s == 0:  # If the vectors are already aligned
        rotation_matrix = np.eye(3)
    else:
        vx = np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]
        ])
        rotation_matrix = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))

    # Apply rotation
    rotated_points = np.dot(translated_points, rotation_matrix.T)

    # Compute inverse rotation matrix
    inverse_rotation = rotation_matrix.T
    return rotated_points, origin, inverse_rotation

def inverse_transform_b1_vector(points, origin, inverse_rotation):
    """Transform points back to the original 3D coordinates."""
    rotated_back_points = np.dot(points, inverse_rotation)
    return rotated_back_points + origin

def calculate_heights_b1_vector(origin, points, hull):
    """Calculate heights from the origin to each edge of the convex hull."""
    edges = hull.simplices  # Edges of the convex hull
    heights = []

    for edge in edges:
        # Points defining the edge
        p1 = points[edge[0]]
        p2 = points[edge[1]]

        # Vector along the edge
        edge_vector = p2 - p1

        # Vector from origin to one point on the edge
        origin_vector = -p1

        # Compute the perpendicular distance from origin to the edge
        height_vector = origin_vector - np.dot(origin_vector, edge_vector) / np.dot(edge_vector, edge_vector) * edge_vector
        height = np.linalg.norm(height_vector)

        # Midpoint of the edge (for visualization)
        midpoint = (p1 + p2) / 2
        heights.append((height, midpoint, p1, p2))

    return heights

def plot_heights_3d(ax, heights, origin, inverse_rotation):
    
    """Plot heights from the origin to each edge in the original 3D coordinates."""
    # Find the shortest height
    shortest_height = min(heights, key=lambda h: h[0])

    for height, midpoint, p1, p2 in heights:
        color = "blue" 
        p1_original = inverse_transform_b1_vector(p1.reshape(1, -1), origin, inverse_rotation)[0]
        p2_original = inverse_transform_b1_vector(p2.reshape(1, -1), origin, inverse_rotation)[0]
        midpoint_original = inverse_transform_b1_vector(midpoint.reshape(1, -1), origin, inverse_rotation)[0]

        ax.plot([p1_original[0], p2_original[0]], [p1_original[1], p2_original[1]], [p1_original[2], p2_original[2]], color="black")
        ax.plot([origin[0], midpoint_original[0]], [origin[1], midpoint_original[1]], [origin[2], midpoint_original[2]], color=color)

def get_b1_vector(atoms, idx1, idx2):

    # Transform the coordinates
    transformed_points, origin, inverse_rotation = transform_coordinates_b1_vector(points, idx1, idx2)

    # Compute convex hull and heights
    hull = ConvexHull(transformed_points[:, :2])
    heights = calculate_heights_b1_vector(np.array([0, 0, 0]), transformed_points, hull)

    shortest_height = min(heights, key=lambda h: h[0])
    height = shortest_height[0]
    midpoint = shortest_height[1]
    midpoint_original = inverse_transform_b1_vector(midpoint.reshape(1, -1), origin, inverse_rotation)[0]
    
    return origin, midpoint_original    
    
    
def on_click(event):
    global selected_points, fig, ax

    idx = find_closest_point(event)
    if idx is not None:
        selected_points.append(idx)
        print(f"Selected point {idx}")

        # Highlight the selected point
        ax.scatter(
            points[idx, 0],
            points[idx, 1],
            points[idx, 2],
            color="green",
            s=100,
        )
        fig.canvas.draw()

        if len(selected_points) == 2:
            idx1, idx2 = selected_points

            # Transform the coordinates
            transformed_points, origin, inverse_rotation = transform_coordinates_b1_vector(points, idx1, idx2)

            # Compute convex hull and heights
            hull = ConvexHull(transformed_points[:, :2])
            heights = calculate_heights_b1_vector(np.array([0, 0, 0]), transformed_points, hull)

            # Replot in original 3D space
            fig.clear()
            ax = fig.add_subplot(111, projection="3d")
            ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="red")

            # Label the points
            for i, (x, y, z) in enumerate(points):
                ax.text(x, y, z, str(i), fontsize=8)

            # Plot transformed edges and heights back in original 3D space
            plot_heights_3d(ax, heights, origin, inverse_rotation)
            
            origin1, midpoint_original1 = get_b1_vector(points, idx1, idx2)

            ax.plot([origin1[0], midpoint_original1[0]], [origin1[1], midpoint_original1[1]], [origin1[2], midpoint_original1[2]], color="Red")

            ax.plot([points[idx1, 0], points[idx2, 0]], [points[idx1, 1], points[idx2, 1]], [points[idx1, 2], points[idx2, 2]], color="purple")

            ax.set_title("Heights and Shape Transformed Back to 3D Space")
            ax.set_xlabel("X-axis")
            ax.set_ylabel("Y-axis")
            ax.set_zlabel("Z-axis")
            plt.show()

# Find the closest point based on mouse click
def find_closest_point(event):
    if event.xdata is None or event.ydata is None:
        return None  # Ignore clicks outside the plot area

    # Map 2D click to 3D space
    x2d, y2d = event.xdata, event.ydata
    distances = []
    for i, (x, y, z) in enumerate(points):
        x3d, y3d, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
        distances.append((i, np.sqrt((x3d - x2d) ** 2 + (y3d - y2d) ** 2)))
    closest_idx = min(distances, key=lambda t: t[1])[0]
    return closest_idx

# Initial plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="red")

# Label each point with its index
for i, (x, y, z) in enumerate(points):
    ax.text(x, y, z, str(i), fontsize=8)

ax.set_title("Click on two points to define new axes")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

selected_points = []  # Store the indices of selected points

# Connect the click event to the callback
from mpl_toolkits.mplot3d import proj3d
fig.canvas.mpl_connect("button_press_event", on_click)

plt.show()
