import numpy as np
from Orbit_conversion import cartesian_to_keplerian

mu = 398600  # Gravitational parameter for Earth [km^3/s^2]

# Define the position and velocity vectors in km and km/s, respectively
r = np.array([10000.0, 40000.0, -5000.0])  # Position vector in km
v = np.array([-1.5, 1.0, -0.1])        # Velocity vector in km/s

# Call the function
keplerian_elements = cartesian_to_keplerian(r, v, mu)

# Print the results
print("Keplerian Elements:")
print(f"Semi-major axis (a): {keplerian_elements[0]} km")
print(f"Eccentricity (e): {keplerian_elements[1]}")
print(f"Inclination (I): {keplerian_elements[2]} degrees")
print(f"Right Ascension of Ascending Node (RAAN): {keplerian_elements[3]} degrees")
print(f"Argument of Periapsis (w): {keplerian_elements[4]} degrees")
print(f"Mean Anomaly (M): {keplerian_elements[5]} degrees")
