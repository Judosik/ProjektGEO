import numpy as np
from scipy.optimize import least_squares

# Original points
points = np.array([(0.139088017115878, -0.117307330471801, -0.545010728916636), (0.0948746468443162, 0.0499956324984577, -0.783813823023455), (0.172694357295053, 0.003001053375046, -0.711388259991645), (-0.0302964750175954, -0.107416583973138, -0.572014800172274), (0.0495682630245263, -0.112120122253242, -0.559431702411444), (-0.108966178147607, -0.0618260937798975, -0.641377338701428), (0.150282872210771, -0.0772939997486353, -0.600264170555548), (0.00328898429940823, 0.0127261026208298, -0.738534121783712), (-0.0974777286959767, -0.0207955873554942, -0.698516112150754), (0.0832030188308752, 0.00816894902498809, -0.725819719509454), (0.161776493899172, -0.0362066878165916, -0.657399193332392), (-0.0190910965034237, -0.0675388012645368, -0.627378251259008), (0.0607873126928254, -0.0721422787648371, -0.61473502840081), (-0.00757433029386827, -0.0263394038640314, -0.684459708475044)])

# Function to fit a plane to points
def fit_plane(points):
    # Reshape points
    points = points.reshape(-1, 3)

    # Create coefficient matrix
    A = np.c_[points[:, 0], points[:, 1], np.ones(points.shape[0])]

    # Solve Ax = b
    normal = np.linalg.lstsq(A, points[:, 2], rcond=None)[0]

    return normal

# Fit plane to points
plane_normal = fit_plane(points)

# Use plane normal as z-axis of new spatial reference system
z_axis = plane_normal / np.linalg.norm(plane_normal)

# Choose arbitrary vectors perpendicular to z-axis as x-axis and y-axis
x_axis = np.array([1, 0, -plane_normal[0] / plane_normal[2]])
x_axis /= np.linalg.norm(x_axis)

y_axis = np.cross(z_axis, x_axis)

# Transformation matrix
transformation_matrix = np.vstack((x_axis, y_axis, z_axis)).T

# Transform points to new spatial reference system
transformed_points = points.dot(transformation_matrix)

print("New spatial reference system:")
print("X-axis:", x_axis)
print("Y-axis:", y_axis)
print("Z-axis:", z_axis)
print("Transformed points:", transformed_points)
