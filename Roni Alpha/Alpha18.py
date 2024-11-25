import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

# Set random seed for reproducibility
np.random.seed(42)

# Generate random 3D points
def generate_points(num_points=20, range_min=-10, range_max=10):
    return pd.DataFrame({
        'x': np.random.uniform(range_min, range_max, num_points),
        'y': np.random.uniform(range_min, range_max, num_points),
        'z': np.random.uniform(range_min, range_max, num_points),
    })

# Select two points by index
def select_points(df):
    print("Available points:\n", df)
    while True:
        try:
            idx_a = int(input("\nSelect index for Point A: "))
            idx_b = int(input("Select index for Point B: "))
            if idx_a == idx_b:
                print("Points A and B must be distinct.")
                continue
            if idx_a not in df.index or idx_b not in df.index:
                print("Invalid indices. Please select valid point indices.")
                continue
            break
        except ValueError:
            print("Please enter valid integer indices.")
    point_a = df.iloc[idx_a].values
    point_b = df.iloc[idx_b].values
    return point_a, point_b

def prepare_data(df):
    """
    Prepares data for finding the shortest vector and visualization.

    Parameters:
        df (DataFrame): DataFrame containing 3D points.

    Returns:
        tuple: Contains point_a, point_b, basis_v, basis_w, projected_points_3d, convex_hull, and ab_vector.
    """
    print("Generated Points:\n", df)

    # Select points A and B
    point_a, point_b = select_points(df)

    # Compute vector AB and basis vectors for the plane
    ab_vector = point_b - point_a
    basis_v, basis_w = calculate_plane_basis(ab_vector)

    # Project points onto the plane
    projected_points = project_to_plane(df.values, point_a, basis_v, basis_w)

    # Compute the convex hull in 2D
    convex_hull = compute_convex_hull(projected_points)

    # Transform projected points back to 3D for visualization
    projected_points_3d = np.array([point_a + p[0] * basis_v + p[1] * basis_w for p in projected_points])

    return point_a, point_b, basis_v, basis_w, projected_points_3d, convex_hull, ab_vector


# Calculate the plane basis vectors perpendicular to AB
def calculate_plane_basis(ab_vector):
    ab_unit = ab_vector / np.linalg.norm(ab_vector)
    arbitrary_vector = np.array([1, 0, 0]) if not np.allclose(ab_unit, [1, 0, 0]) else np.array([0, 1, 0])
    basis_v = np.cross(ab_unit, arbitrary_vector)
    basis_v /= np.linalg.norm(basis_v)
    basis_w = np.cross(ab_unit, basis_v)
    return basis_v, basis_w

# Project points onto the plane
def project_to_plane(points, point_a, basis_v, basis_w):
    projected_points = []
    for point in points:
        ap = point - point_a
        x = np.dot(ap, basis_v)
        y = np.dot(ap, basis_w)
        projected_points.append([x, y])
    return np.array(projected_points)

# Find the convex hull
def compute_convex_hull(points_2d):
    return ConvexHull(points_2d)

# Find the shortest distance from a point to a segment
def distance_to_segment(point, segment_start, segment_end):
    segment_vector = segment_end - segment_start
    segment_length = np.linalg.norm(segment_vector)
    if segment_length == 0:
        return np.linalg.norm(point - segment_start), segment_start
    t = np.dot(point - segment_start, segment_vector) / segment_length**2
    t = np.clip(t, 0, 1)
    projection = segment_start + t * segment_vector
    distance = np.linalg.norm(point - projection)
    return distance, projection

# Find the shortest vector from A to any edge of the convex hull in 2D
def find_shortest_vector_2d(point_a_2d, convex_hull):
    min_distance = float('inf')
    closest_projection = None
    for i in range(len(convex_hull.vertices)):
        start = convex_hull.points[convex_hull.vertices[i]]
        end = convex_hull.points[convex_hull.vertices[(i + 1) % len(convex_hull.vertices)]]
        distance, projection = distance_to_segment(point_a_2d, start, end)
        if distance < min_distance:
            min_distance = distance
            closest_projection = projection
    return closest_projection

# Compute the shortest vector in 3D
def get_shortest_vector_3d(df, point_a, point_b):
    ab_vector = point_b - point_a
    basis_v, basis_w = calculate_plane_basis(ab_vector)

    # Project points onto the plane
    projected_points = project_to_plane(df.values, point_a, basis_v, basis_w)

    # Compute convex hull in 2D
    convex_hull = compute_convex_hull(projected_points)

    # Find the shortest vector in 2D
    point_a_2d = np.array([0, 0])  # Point A in the plane
    closest_proj_2d = find_shortest_vector_2d(point_a_2d, convex_hull)

    # Back-project the shortest vector to 3D
    if closest_proj_2d is not None:
        shortest_vector_3d = closest_proj_2d[0] * basis_v + closest_proj_2d[1] * basis_w
    else:
        shortest_vector_3d = np.array([0, 0, 0])  # Default to zero vector if no projection found
    print(shortest_vector_3d)
    print("shortest_vector_3d")
    return shortest_vector_3d
   

# Visualize in 3D
def visualize_3d(df, point_a, point_b, projected_points_3d, convex_hull, shortest_vector_3d):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_proj_type('ortho')

    # Plot original points
    ax.scatter(df['x'], df['y'], df['z'], label='Original Points', alpha=0.6)

    # Plot points A and B
    ax.scatter(*point_a, color='green', label='Point A', s=100)
    ax.scatter(*point_b, color='red', label='Point B', s=100)

    # Plot projected points on the plane
    ax.scatter(projected_points_3d[:, 0], projected_points_3d[:, 1], projected_points_3d[:, 2], label='Projected Points', alpha=0.8)

    # Plot convex hull edges
    for simplex in convex_hull.simplices:
        start = projected_points_3d[simplex[0]]
        end = projected_points_3d[simplex[1]]
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color='purple', label='Convex Hull Edge' if simplex[0] == 0 else "")

    # Plot shortest vector
    ax.quiver(point_a[0], point_a[1], point_a[2],
              shortest_vector_3d[0], shortest_vector_3d[1], shortest_vector_3d[2],
              color='orange', label='Shortest Vector', linewidth=2)

    ax.legend()
    ax.set_title("3D Visualization with Orthographic Camera")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    plt.show()

# Main program workflow
df = generate_points()
print("Generated Points:\n", df)

point_a, point_b = select_points(df)

# Compute the shortest vector and related data for visualization
# Compute the shortest vector and related data for visualization
ab_vector = point_b - point_a
basis_v, basis_w = calculate_plane_basis(ab_vector)
projected_points = project_to_plane(df.values, point_a, basis_v, basis_w)
convex_hull = compute_convex_hull(projected_points)
shortest_vector_3d = get_shortest_vector_3d(df, point_a, point_b)

# Transform 2D projected points back to 3D
projected_points_3d = np.array([point_a + p[0] * basis_v + p[1] * basis_w for p in projected_points])

# Visualize results
visualize_3d(df, point_a, point_b, projected_points_3d, convex_hull, shortest_vector_3d)


# Visualize results
visualize_3d(df, point_a, point_b, np.array(projected_points), convex_hull, shortest_vector_3d)
