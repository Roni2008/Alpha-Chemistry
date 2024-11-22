def calculate_plane_basis(ab_vector):
    ab_unit = ab_vector / np.linalg.norm(ab_vector)
    arbitrary_vector = np.array([1, 0, 0]) if not np.allclose(ab_unit, [1, 0, 0]) else np.array([0, 1, 0])
    basis_v = np.cross(ab_unit, arbitrary_vector)
    basis_v /= np.linalg.norm(basis_v)
    basis_w = np.cross(ab_unit, basis_v)
    return basis_v, basis_w

# Project points onto the plane
def project_to_plane(points, point_a, ab_vector, basis_v, basis_w):
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

# Main pre-computation workflow
df = generate_points()
print("Generated Points:\n", df)

point_a, point_b = select_points(df)
ab_vector = point_b - point_a
basis_v, basis_w = calculate_plane_basis(ab_vector)

# Project points onto the plane
projected_points = project_to_plane(df.values, point_a, ab_vector, basis_v, basis_w)

# Compute convex hull in 2D
convex_hull = compute_convex_hull(projected_points)

# Find the shortest vector in 2D
point_a_2d = np.array([0, 0])
shortest_vector_2d = None  # Placeholder for global use later

# Back-project 2D points and vectors to 3D
projected_points_3d = [point_a + p[0] * basis_v + p[1] * basis_w for p in projected_points]
shortest_vector_3d = None  # Placeholder for global use later

# Find the shortest vector from A to any edge of the convex hull
def get_shortest_vector(point_a_2d, convex_hull):
    """
    Finds the shortest vector from a point to any edge of the convex hull.

    Parameters:
        point_a_2d (np.array): The 2D coordinates of the reference point.
        convex_hull (ConvexHull): The convex hull object computed from points.

    Returns:
        np.array: The shortest vector from the point to the closest edge of the convex hull.
    """
    min_distance = float('inf')
    closest_projection = None
    for i in range(len(convex_hull.vertices)):
        start = convex_hull.points[convex_hull.vertices[i]]
        end = convex_hull.points[convex_hull.vertices[(i + 1) % len(convex_hull.vertices)]]
        distance, projection = distance_to_segment(point_a_2d, start, end)
        if distance < min_distance:
            min_distance = distance
            closest_projection = projection
    shortest_vector = closest_projection - point_a_2d
    return shortest_vector
