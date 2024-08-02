import numpy as np
from scipy.optimize import newton
import csv

# Constants
G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
M_earth = 5.972e24  # Mass of Earth in kg
mu = G * M_earth  # Gravitational parameter for Earth in m^3/s^2

# Provided functions
def mean_to_eccentric_anomaly(M, e, tol=1e-10):
    """Solve for Eccentric Anomaly (E) given Mean Anomaly (M) and eccentricity (e)"""
    func = lambda E: E - e * np.sin(E) - M
    E0 = M if e < 0.8 else np.pi
    return newton(func, E0, tol=tol)

def keplerian_to_cartesian(a_km, e, I, Omega, omega, M, mu):
    """Convert from Keplerian elements to Cartesian state vector."""
    # Convert semi-major axis to meters for calculation
    a = a_km * 1000  # Convert km to m

    I, Omega, omega = np.radians(I), np.radians(Omega), np.radians(omega)

    E = mean_to_eccentric_anomaly(M, e)

    nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

    r_peri = a * (1 - e * np.cos(E))
    r_vec_peri = np.array([r_peri * np.cos(nu), r_peri * np.sin(nu), 0])

    v_vec_peri = np.array([-np.sin(E), np.sqrt(1 - e**2) * np.cos(E), 0])
    v_vec_peri *= np.sqrt(mu * a) / r_peri / np.sqrt(1 - e**2)

    R_z_Omega = np.array([[np.cos(Omega), -np.sin(Omega), 0], [np.sin(Omega), np.cos(Omega), 0], [0, 0, 1]])
    R_x_I = np.array([[1, 0, 0], [0, np.cos(I), -np.sin(I)], [0, np.sin(I), np.cos(I)]])
    R_z_omega = np.array([[np.cos(omega), -np.sin(omega), 0], [np.sin(omega), np.cos(omega), 0], [0, 0, 1]])
    R = R_z_Omega @ R_x_I @ R_z_omega

    r_vec = R @ r_vec_peri
    v_vec = R @ v_vec_peri

    # Convert results to km and km/s for output
    return r_vec / 1000, v_vec / 1000

# Function to calculate the orbital period
def orbital_period(a_km):
    a = a_km * 1000  # Convert km to m for calculation
    return 2 * np.pi * np.sqrt(a**3 / mu)

# Function to predict the orbit for 16 revolutions
def predict_orbit_and_write_csv(a_km, e, I, Omega, omega, nu, filename="satellite_orbit.csv", num_points_per_orbit=100):
    T = orbital_period(a_km)  # Orbital period
    delta_t = T / num_points_per_orbit  # Time step

    orbit_positions = []

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["X (km)", "Y (km)", "Z (km)", "Time Step", "Revolution"])

        for revolution in range(16):
            for step in range(num_points_per_orbit):
                M = 2 * np.pi * step * delta_t / T  # Mean anomaly
                r_vec, _ = keplerian_to_cartesian(a_km, e, I, Omega, omega, M, mu)
                orbit_positions.append(r_vec)

                # Write to CSV: position (in km) and time information
                writer.writerow([r_vec[0], r_vec[1], r_vec[2], step, revolution])

    return np.array(orbit_positions)



if __name__ == "__main__":
    # Example orbital elements (replace with actual values)
    a = 6827.9986  # Semi-major axis in kilometers (including Earth's radius if altitude is given)
    e = 0.004    # Eccentricity
    I = 87.3    # Inclination in degrees
    Omega = 0  # Right Ascension of Ascending Node in degrees
    omega = 0  # Argument of Perigee in degrees
    nu = 0.2106   # True Anomaly in degrees

    # Predict the orbit and write to CSV
    orbit_positions = predict_orbit_and_write_csv(a, e, I, Omega, omega, nu, r"C:\Users\KhushaldasBadhan\OneDrive - ODYSSEUS SPACE\Documents\semester\GNSS\HW3\predicted_orbit.csv")

    # Show the shape of the orbit positions array
    orbit_positions.shape  # This should give us 16 orbits with 100 points each

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Function to plot the orbit
    def plot_orbit(orbit_positions):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the orbit
        ax.plot(orbit_positions[:, 0], orbit_positions[:, 1], orbit_positions[:, 2], label='Satellite Orbit')

        # Plot Earth (assuming a perfect sphere for simplicity)
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = Earth_radius * np.cos(u) * np.sin(v)
        y = Earth_radius * np.sin(u) * np.sin(v)
        z = Earth_radius * np.cos(v)
        ax.plot_surface(x, y, z, color='b', alpha=0.3)

        # Axes labels
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('3D Orbit Plot')

        plt.legend()
        plt.show()

    # Radius of the Earth
    Earth_radius = 6371e3  # in meters

    # Plot the orbit
    plot_orbit(orbit_positions)

