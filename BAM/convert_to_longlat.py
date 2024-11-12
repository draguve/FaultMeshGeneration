from pyproj import Transformer
import math

# Define transformers for each projection system
transformer_utm = Transformer.from_crs("EPSG:32610", "EPSG:4326", always_xy=True)  # UTM Zone 10N (WGS 84)

# Coordinates (Easting, Northing) format for each point
points = [
    (437974.13, 4284863.14),  # Point 1
    (706265.49, 4284863.14),  # Point 2
    (706265.49, 4029171.34),  # Point 3
    (437974.13, 4029171.34)   # Point 4
]

# Transform points using UTM Zone 10N
transformed_points_utm = [transformer_utm.transform(*point) for point in points]
# Format as MULTILINESTRING for UTM Zone 10N
utm_multilinestring = "MULTILINESTRING((" + ", ".join(f"{lon} {lat} 0.0" for lon, lat in transformed_points_utm) + "))"

# Output the formatted strings
print("Using UTM Zone 10N (WGS 84) - EPSG:32610:")
print(utm_multilinestring)

# Define transformers for each projection system
transformer_utm = Transformer.from_crs("EPSG:32610", "EPSG:4326", always_xy=True)  # UTM Zone 10N (WGS 84)

# Coordinates (Easting, Northing) format for each point
points = [
    (437974.13, 4284863.14),  # Point 1
    (706265.49, 4284863.14),  # Point 2
    (706265.49, 4029171.34),  # Point 3
    (437974.13, 4029171.34)   # Point 4
]

# Transform points using UTM Zone 10N
transformed_points_utm = [transformer_utm.transform(*point) for point in points]

# Rotation angle in degrees
theta_deg = 35+90
theta_rad = math.radians(theta_deg)

# Rotation matrix components
cos_theta = math.cos(theta_rad)
sin_theta = math.sin(theta_rad)

# Function to rotate points
def rotate_point(lon, lat, origin_lon, origin_lat):
    # Translate the point to origin
    x = lon - origin_lon
    y = lat - origin_lat
    # Apply the rotation matrix
    x_rot = cos_theta * x - sin_theta * y
    y_rot = sin_theta * x + cos_theta * y
    # Translate back to original position
    lon_rot = x_rot + origin_lon
    lat_rot = y_rot + origin_lat
    return lon_rot, lat_rot

# Get the centroid (average of the points) for rotation origin
avg_lon = sum(lon for lon, lat in transformed_points_utm) / len(transformed_points_utm)
avg_lat = sum(lat for lon, lat in transformed_points_utm) / len(transformed_points_utm)

# Rotate each point
rotated_points = [rotate_point(lon, lat, avg_lon, avg_lat) for lon, lat in transformed_points_utm]

# Format as MULTILINESTRING for rotated points
utm_multilinestring = "MULTILINESTRING((" + ", ".join(f"{lon} {lat} 0.0" for lon, lat in rotated_points) + "))"

# Output the formatted strings
print("Using UTM Zone 10N (WGS 84) - EPSG:32610:")
print(utm_multilinestring)


# import matplotlib.pyplot as plt
# from pyproj import Transformer
#
# # Define transformers for each projection system
# transformer_utm = Transformer.from_crs("EPSG:32610", "EPSG:4326", always_xy=True)  # UTM Zone 10N (WGS 84)
# transformer_ca = Transformer.from_crs("EPSG:2227", "EPSG:4326", always_xy=True)    # NAD83 / California Zone 3 (ftUS)
#
# # Coordinates (Easting, Northing) format for each point
# points = [
#     (437974.13, 4284863.14),  # Point 1
#     (706265.49, 4284863.14),  # Point 2
#     (706265.49, 4029171.34),  # Point 3
#     (437974.13, 4029171.34)   # Point 4
# ]
#
# # Transform points using UTM Zone 10N
# transformed_points_utm = [transformer_utm.transform(*point) for point in points]
# # Transform points using NAD83 / California Zone 3
# transformed_points_ca = [transformer_ca.transform(*point) for point in points]
#
# # Plotting the points
# plt.figure(figsize=(10, 5))
#
# # Plot UTM transformed points
# utm_lons, utm_lats = zip(*transformed_points_utm)
# plt.plot(utm_lons, utm_lats, 'o-', color='blue', label="UTM Zone 10N (WGS 84)")
# for i, (lon, lat) in enumerate(zip(utm_lons, utm_lats), 1):
#     plt.text(lon, lat, f"Point {i}", fontsize=9, color='blue')
#
# # # Plot California Zone 3 transformed points
# # ca_lons, ca_lats = zip(*transformed_points_ca)
# # plt.plot(ca_lons, ca_lats, 'o-', color='red', label="NAD83 / California Zone 3")
# # for i, (lon, lat) in enumerate(zip(ca_lons, ca_lats), 1):
# #     plt.text(lon, lat, f"Point {i}", fontsize=9, color='red')
#
# # Adding labels and legend
# plt.xlabel("Longitude")
# plt.ylabel("Latitude")
# plt.legend()
# plt.title("Transformed Points in Different Coordinate Systems")
# plt.grid(True)
# plt.show()