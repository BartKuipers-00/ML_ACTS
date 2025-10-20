import numpy as np
import matplotlib.pyplot as plt
import math

def bethe_bloch_dEdx(p_GeV, mass_GeV=0.140, Z=1, A=1, rho_g_cm3=1.0, I_eV=75):
    """
    Simplified Bethe-Bloch energy loss formula
    
    Parameters:
    - p_GeV: momentum in GeV/c
    - mass_GeV: particle mass in GeV/c² (default: pion = 0.140)
    - Z: atomic number of material (default: 1 for hydrogen)
    - A: atomic mass of material (default: 1)
    - rho_g_cm3: material density in g/cm³
    - I_eV: mean excitation energy in eV (default: ~75 eV for light materials)
    
    Returns: dE/dx in MeV/cm
    """
    
    # Constants
    K = 0.307075  # MeV·cm²/mol for A=1
    me = 0.000511  # electron mass in GeV
    
    # Calculate gamma and beta
    E = np.sqrt(p_GeV**2 + mass_GeV**2)  # Total energy
    gamma = E / mass_GeV
    beta = p_GeV / E
    
    # Avoid numerical issues
    if beta < 0.01:
        return 1000.0  # Very high energy loss for very slow particles
    
    # Mean excitation energy in GeV
    I_GeV = I_eV * 1e-9
    
    # Bethe-Bloch formula (simplified, no density correction)
    Wmax = 2 * me * beta**2 * gamma**2 / (1 + 2*gamma*me/mass_GeV + (me/mass_GeV)**2)
    
    ln_term = np.log(2 * me * beta**2 * gamma**2 * Wmax / I_GeV**2)
    
    dEdx = K * (Z/A) * rho_g_cm3 * (1/beta**2) * (0.5 * ln_term - beta**2)
    
    return max(dEdx, 0.0)  # Ensure positive energy loss


def p_to_pT(p_GeV, theta_deg):
    """Convert total momentum to transverse momentum"""
    theta_rad = math.radians(float(theta_deg))
    return p_GeV * math.sin(theta_rad)


def radius_from_pT(pT_GeV, B_T=3.8, q=-1.0):
    """Calculate radius of curvature from transverse momentum"""
    if pT_GeV <= 0:
        return 0.0
    return pT_GeV / (0.3 * abs(q) * B_T)


def trajectory_with_energy_loss(p0_GeV, theta_deg, B_T=3.8, q=-1.0, 
                               mass_GeV=0.140, material_params=None,
                               max_distance=1.0, step_size_cm=0.5):
    """
    Calculate particle trajectory accounting for energy loss
    
    Parameters:
    - p0_GeV: initial momentum in GeV/c
    - theta_deg: polar angle in degrees
    - B_T: magnetic field in Tesla
    - q: charge in elementary charges
    - mass_GeV: particle mass in GeV/c²
    - material_params: dict with Z, A, rho_g_cm3, I_eV for material
    - max_distance: maximum distance to calculate in meters
    - step_size_cm: step size for integration in cm
    
    Returns: (distances_m, momenta_GeV, radii_m, stopped_flag)
    """
    
    # Default material parameters (silicon-like)
    if material_params is None:
        material_params = {
            'Z': 14,      # Silicon
            'A': 28,
            'rho_g_cm3': 2.33,
            'I_eV': 173
        }
    
    # Initialize arrays
    distances = [0.0]  # cm
    momenta = [p0_GeV]
    radii = []
    
    current_p = p0_GeV
    current_distance = 0.0
    
    while current_distance < max_distance * 100:  # Convert to cm
        # Calculate current transverse momentum and radius
        pT = p_to_pT(current_p, theta_deg)
        radius = radius_from_pT(pT, B_T, q)
        radii.append(radius)
        
        # Check if particle stopped (momentum too low)
        if current_p < 0.01:  # 10 MeV threshold
            break
            
        # Calculate energy loss rate
        dEdx = bethe_bloch_dEdx(current_p, mass_GeV, **material_params)
        
        # Convert dE/dx to dp/dx (approximately dE ≈ dp for relativistic particles)
        # More precisely: dE = p*dp/sqrt(p² + m²) ≈ dp for p >> m
        E = np.sqrt(current_p**2 + mass_GeV**2)
        dpdx = dEdx * current_p / E  # GeV/cm
        
        # Update momentum
        dp = dpdx * step_size_cm
        current_p = max(current_p - dp, 0.0)
        
        # Update distance
        current_distance += step_size_cm
        
        distances.append(current_distance)
        momenta.append(current_p)
    
    # Convert distances to meters
    distances_m = np.array(distances) / 100.0
    momenta_GeV = np.array(momenta)
    radii_m = np.array(radii)
    
    stopped_flag = current_p < 0.01
    
    return distances_m, momenta_GeV, radii_m, stopped_flag


def max_radial_reach_with_energy_loss(p0_GeV, theta_deg, **kwargs):
    """
    Calculate maximum radial distance a particle reaches before stopping
    or completing a semicircle, accounting for energy loss
    """
    distances, momenta, radii, stopped = trajectory_with_energy_loss(
        p0_GeV, theta_deg, **kwargs
    )
    
    if len(radii) == 0:
        return 0.0
        
    # Maximum radial reach is approximately the diameter at the point
    # where the particle either stops or completes its trajectory
    if stopped:
        # Particle stopped - use final radius * 2
        return 2 * radii[-1] if radii[-1] > 0 else 0.0
    else:
        # Particle didn't stop in our calculation range
        return 2 * radii[0]  # Initial diameter


# Test the implementation
if __name__ == "__main__":
    # Test single trajectory
    p0 = 0.2  # GeV
    theta = 90  # degrees (eta = 0)
    
    distances, momenta, radii, stopped = trajectory_with_energy_loss(
        p0, theta, B_T=3.8, max_distance=0.5
    )
    
    print(f"Initial momentum: {p0:.3f} GeV")
    print(f"Theta: {theta}°")
    print(f"Stopped: {stopped}")
    print(f"Final momentum: {momenta[-1]:.3f} GeV")
    print(f"Distance traveled: {distances[-1]:.3f} m")
    print(f"Initial radius: {radii[0]:.3f} m")
    if len(radii) > 1:
        print(f"Final radius: {radii[-1]:.3f} m")
    
    # Plot trajectory evolution
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Momentum vs distance
    ax1.plot(distances, momenta)
    ax1.set_xlabel('Distance [m]')
    ax1.set_ylabel('Momentum [GeV/c]')
    ax1.set_title('Momentum vs Distance')
    ax1.grid(True)
    
    # Radius vs distance
    ax2.plot(distances[:-1], radii)  # radii has one less element
    ax2.set_xlabel('Distance [m]')
    ax2.set_ylabel('Radius [m]')
    ax2.set_title('Radius of Curvature vs Distance')
    ax2.grid(True)
    
    # Energy loss rate vs momentum
    p_range = np.linspace(0.05, 1.0, 100)
    dedx_values = [bethe_bloch_dEdx(p) for p in p_range]
    ax3.semilogy(p_range, dedx_values)
    ax3.set_xlabel('Momentum [GeV/c]')
    ax3.set_ylabel('dE/dx [MeV/cm]')
    ax3.set_title('Bethe-Bloch Energy Loss')
    ax3.grid(True)
    
    # Comparison: with vs without energy loss
    theta_test = 90
    p_range_test = np.linspace(0.05, 0.5, 50)
    
    # Without energy loss (original formula)
    reach_no_loss = [2 * p_to_pT(p, theta_test) / (0.3 * 3.8) for p in p_range_test]
    
    # With energy loss
    reach_with_loss = [max_radial_reach_with_energy_loss(p, theta_test, B_T=3.8) 
                       for p in p_range_test]
    
    ax4.plot(p_range_test, reach_no_loss, 'b-', label='No energy loss')
    ax4.plot(p_range_test, reach_with_loss, 'r-', label='With energy loss')
    ax4.axhline(y=0.26, color='k', linestyle='--', label='Strip radius (0.26 m)')
    ax4.set_xlabel('Initial Momentum [GeV/c]')
    ax4.set_ylabel('Maximum Radial Reach [m]')
    ax4.set_title(f'Radial Reach vs Momentum (θ={theta_test}°)')
    ax4.legend()
    ax4.grid(True)
    
    plt.tight_layout()
    plt.show()