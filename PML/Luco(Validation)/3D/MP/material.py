# %%
def calculate_elastic_and_shear_modulus(density, poisson_ratio, shear_wave_velocity):
    # Calculate shear modulus
    shear_modulus = density * shear_wave_velocity**2 

    # Calculate elastic modulus (Young's modulus)
    elastic_modulus = 2 * shear_modulus * (1 + poisson_ratio)
    

    return elastic_modulus, shear_modulus

# Example usage:
density = 1800  # kg/m^3 (density)
poisson_ratio = 0.25  # (Poisson's ratio)
shear_wave_velocity = 300  # m/s (shear wave velocity)

elastic_modulus, shear_modulus = calculate_elastic_and_shear_modulus(density, poisson_ratio, shear_wave_velocity)
print(f"Elastic Modulus: {elastic_modulus} Pa")
print(f"Shear Modulus: {shear_modulus} Pa")


# %%
