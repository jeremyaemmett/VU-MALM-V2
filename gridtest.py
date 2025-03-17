import numpy as np
import matplotlib.pyplot as plt


def convert_to_exponential_grid(linear_grid, gas_concentrations, num_exp_points):
    """
    Converts gas concentrations on a high-resolution linear spaced grid to a lower-resolution exponentially spaced grid
    while conserving mass and shape.

    Parameters:
    - linear_grid: A 1D numpy array representing the high-resolution linear grid (e.g., distances, time, etc.).
    - gas_concentrations: A 1D numpy array representing the gas concentrations corresponding to the linear grid.
    - num_exp_points: The number of points for the exponentially spaced grid.

    Returns:
    - exp_grid: The exponentially spaced grid.
    - exp_concentrations: The gas concentrations on the exponentially spaced grid.
    """
    # Step 1: Create an exponentially spaced grid (ranging from start to end of the linear grid)
    start = linear_grid[0]
    end = linear_grid[-1]
    exp_grid = np.logspace(np.log10(start), np.log10(end), num_exp_points)

    # Step 2: Interpolate the gas concentrations onto the exponentially spaced grid
    exp_concentrations = np.interp(exp_grid, linear_grid, gas_concentrations)

    # Step 3: Mass conservation by scaling the concentrations
    linear_mass = np.trapz(gas_concentrations, linear_grid)  # Integral over the linear grid
    exp_mass = np.trapz(exp_concentrations, exp_grid)  # Integral over the exponential grid

    # Step 4: Scale the concentrations on the exponential grid to conserve mass
    scaling_factor = linear_mass / exp_mass
    exp_concentrations *= scaling_factor

    return exp_grid, exp_concentrations


# Example usage
linear_grid = np.linspace(0.1, 10, 1000)  # High-resolution linear grid
gas_concentrations = np.exp(-linear_grid)  # Example gas concentration profile
num_exp_points = 50  # Lower resolution with 50 points

exp_grid, exp_concentrations = convert_to_exponential_grid(linear_grid, gas_concentrations, num_exp_points)

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(linear_grid, gas_concentrations, label='Original Linear Grid', color='blue')
plt.plot(exp_grid, exp_concentrations, label='Converted Exponential Grid', color='red', marker='o')
plt.xscale('log')  # Use logarithmic scale for x-axis to show exponential spacing
plt.xlabel('Grid points')
plt.ylabel('Gas Concentration')
plt.legend()
plt.title('Gas Concentrations on Linear vs Exponential Grid')
plt.show()
