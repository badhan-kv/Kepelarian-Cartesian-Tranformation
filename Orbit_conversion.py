import numpy as np
from scipy.optimize import newton

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

def cartesian_to_keplerian(r, v, mu):
    # Specific angular momentum
    h = np.cross(r, v)
    n = np.cross([0, 0, 1], h)

    # Eccentricity vector
    e_vector = np.cross(v, h) / mu - r / np.linalg.norm(r)
    e = np.linalg.norm(e_vector)

    # Semi-major axis
    energy = np.linalg.norm(v)**2 / 2 - mu / np.linalg.norm(r)
    a = -mu / (2 * energy)

    # Inclination
    I = np.degrees(np.arccos(h[2] / np.linalg.norm(h)))  # Converted to degrees

    # RAAN
    RAAN = np.degrees(np.arccos(n[0] / np.linalg.norm(n)))  # Converted to degrees
    if n[1] < 0:
        RAAN = 360 - RAAN

    # Argument of periapsis
    w = np.degrees(np.arccos(np.dot(n, e_vector) / (np.linalg.norm(n) * e)))  # Converted to degrees
    if e_vector[2] < 0:
        w = 360 - w

    # True anomaly
    nu = np.degrees(np.arccos(np.dot(e_vector, r) / (e * np.linalg.norm(r))))  # Converted to degrees
    if np.dot(r, v) < 0:
        nu = 360 - nu
    
    # Mean Anomaly (M)
    E = 2 * np.arctan(np.tan(nu / 2) / np.sqrt((1 + e) / (1 - e)))
    M = np.degrees(E - e * np.sin(E))  # Converted to degrees
    
    # Normalize M to be within 0 to 360 degrees
    M = M % 360


    return a, e, I, RAAN, w, M

if __name__ == "__main__":
    # Example usage
    mu = 398600  # Gravitational parameter for Earth [km^3/s^2]
    a, e, I, Omega, omega, M = 7000, 0.001, 30, 40, 50, 0  # Keplerian elements

    r_vec, v_vec = keplerian_to_cartesian(a, e, I, Omega, omega, M, mu)
    print("Cartesian State Vector:")
    print("Position Vector (r):", r_vec)
    print("Velocity Vector (v):", v_vec)

    # Convert back from Cartesian to Keplerian
    a_back, e_back, I_back, Omega_back, omega_back, M_back = cartesian_to_keplerian(r_vec, v_vec, mu)
    print("\nConverted Back to Keplerian Elements:")
    print("Semi-major axis (a):", a_back)
    print("Eccentricity (e):", e_back)
    print("Inclination (I):", I_back)
    print("RAAN (Omega):", Omega_back)
    print("Argument of Periapsis (omega):", omega_back)
    print("Mean Anomaly (M):", M_back)