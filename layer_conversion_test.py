import numpy as np
from scipy.integrate import trapz
from scipy.interpolate import interp1d

def convert_to_exponential_spacing(original_concentrations, original_spacing, growth_rate, n_layers):
    # Create the original position array (uniformly spaced)
    original_positions = np.arange(0.05, len(original_concentrations) * original_spacing, original_spacing)

    # Create the new position array (exponentially spaced)
    new_positions = np.array([original_positions[0] * np.exp(growth_rate * i) for i in range(n_layers)])

    # Interpolate the original concentrations to the new positions
    interpolator = interp1d(original_positions, original_concentrations, kind='linear', fill_value='extrapolate')
    new_concentrations = interpolator(new_positions)

    # Compute the total mass of the original concentration profile
    original_mass = trapz(original_concentrations, original_positions)

    # Compute the new concentration array such that mass is conserved
    # First, compute the area under the curve for the exponentially spaced layers
    new_mass = trapz(new_concentrations, new_positions)

    # Scale the concentrations so that the total mass is conserved
    scaling_factor = original_mass / new_mass
    new_concentrations *= scaling_factor

    return new_concentrations, new_positions

# Example usage
original_concentrations = np.array([10, 20, 30, 40, 50])  # example concentrations
original_spacing = 1.0  # uniform spacing between layers
growth_rate = 0.1  # exponential growth rate for spacing
n_layers = 10  # number of layers for the new exponential spacing

# Convert to new concentrations with exponential spacing
new_concentrations, new_positions = convert_to_exponential_spacing(original_concentrations, original_spacing, growth_rate, n_layers)

print("Original concentrations:", original_concentrations)
print("New concentrations with exponential spacing:", new_concentrations)
print("New layer positions:", new_positions)
