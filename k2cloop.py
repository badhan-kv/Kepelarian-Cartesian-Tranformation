import numpy as np
from scipy.optimize import newton

# Function Definitions
def mean_to_eccentric_anomaly(M, e, tol=1e-10):
    """Solve for Eccentric Anomaly (E) given Mean Anomaly (M) and eccentricity (e)"""
    func = lambda E: E - e * np.sin(E) - M
    E0 = M if e < 0.8 else np.pi
    return newton(func, E0, tol=tol)

def keplerian_to_cartesian(a, e, I, Omega, omega, M, mu):
    """Convert from Keplerian elements to Cartesian state vector."""
    # Convert to radians
    I, Omega, omega = np.radians(I), np.radians(Omega), np.radians(omega)

    # Calculate the Eccentric Anomaly
    E = mean_to_eccentric_anomaly(M, e)

    # True Anomaly
    nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

    # Perifocal coordinates
    r_peri = a * (1 - e * np.cos(E))
    r_vec_peri = np.array([r_peri * np.cos(nu), r_peri * np.sin(nu), 0])

    # Adjusting velocity vector calculation
    v_vec_peri = np.array([-np.sin(E), np.sqrt(1 - e**2) * np.cos(E), 0])
    v_vec_peri *= np.sqrt(mu * a) / r_peri / np.sqrt(1 - e**2)

    # Rotation Matrices
    R_z_Omega = np.array([[np.cos(Omega), -np.sin(Omega), 0], [np.sin(Omega), np.cos(Omega), 0], [0, 0, 1]])
    R_x_I = np.array([[1, 0, 0], [0, np.cos(I), -np.sin(I)], [0, np.sin(I), np.cos(I)]])
    R_z_omega = np.array([[np.cos(omega), -np.sin(omega), 0], [np.sin(omega), np.cos(omega), 0], [0, 0, 1]])
    R = R_z_Omega @ R_x_I @ R_z_omega

    # Transform to Ecliptic Coordinates
    r_vec = R @ r_vec_peri
    v_vec = R @ v_vec_peri

    return r_vec, v_vec


# Gravitational parameter for Earth [km^3/s^2] and other orbital parameters
mu = 398600  
a = 10000    # Semi-major axis [km]
e = 0.85     # Eccentricity
I = 0        # Inclination [degrees]
Omega = 45   # Right Ascension of Ascending Node [degrees]
omega = 30   # Argument of Perigee [degrees]

# Prepare to store orbit points
M_values = np.linspace(0, 2*np.pi, 360)  # 360 points from 0 to 2π radians
orbit_points = np.zeros((len(M_values), 3))  # To store orbit points

# Calculate position for each M value
for idx, M in enumerate(M_values):
    r_vec, _ = keplerian_to_cartesian(a, e, I, Omega, omega, M, mu)
    orbit_points[idx] = r_vec

# Save to CSV file
np.savetxt("orbit_data.csv", orbit_points, delimiter=",")

print("done!")
