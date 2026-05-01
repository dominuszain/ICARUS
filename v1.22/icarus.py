#!/usr/bin/env python3
"""
ICARUS v1.22 - Nuclear Reaction Rate Calculator
Python implementation by Zain Ul Abideen
Original Fortran code upgraded to Python with numpy/scipy
"""

import argparse
import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import sys

# Physical constants
PI = 3.14159265358979
C = 2.997e8  # speed of light (m/s)
AVO = 6.022e23  # Avogadro's number (mol^-1)
MASS_NEUTRON = 939.6  # MeV/c^2
AMU = 931.5  # MeV/c^2


def print_logo():
    """Print ICARUS ASCII logo"""
    print("_________ _______  _______  _______           _______ ")
    print("\\__   __/(  ____ \\(  ___  )(  ____ )|\\     /|(  ____ \\")
    print("   ) (   | (    \\/| (   ) || (    )|| )   ( || (    \\/")
    print("   | |   | |      | (___) || (____)|| |   | || (_____ ")
    print("   | |   | |      |  ___  ||     __)| |   | |(_____  )")
    print("   | |   | |      | (   ) || (\\ (   | |   | |      ) |")
    print("___) (___| (____/\\| )   ( || ) \\ \\__| (___) |/\\____) |")
    print("\\_______/(_______/|/     \\||/   \\__/(_______)\\_______) v1.22")
    print("")
    print("Authored by Zain Ul Abideen as part of FYP (BS-SS09-IST)")
    print("contact: dominuszain@gmail.com")
    print("")


def load_cross_section_data(filename):
    """
    Load energy and cross-section data from file.
    
    Parameters
    ----------
    filename : str
        Path to data file
    
    Returns
    -------
    energy : ndarray
        Energy values in MeV
    cross_section : ndarray
        Cross-section values in millibarns
    n_points : int
        Number of data points loaded
    """
    data = np.loadtxt(filename)
    energy = data[:, 0]  # MeV
    cross_section = data[:, 1]  # millibarns
    n_points = len(energy)
    return energy, cross_section, n_points


def calculate_reduced_mass(target_mass):
    """
    Calculate reduced mass for neutron-target system.
    
    Parameters
    ----------
    target_mass : float
        Target nucleus mass number (A)
    
    Returns
    -------
    float
        Reduced mass in MeV/c^2
    """
    return (MASS_NEUTRON * target_mass * AMU) / (target_mass * AMU + MASS_NEUTRON)


def calculate_reaction_rate(incident_energy_kev, energy_data, xs_data, reduced_mass):
    """
    Calculate Maxwellian-averaged cross-section and reaction rate.
    
    Parameters
    ----------
    incident_energy_kev : float
        Incident neutron energy in keV
    energy_data : ndarray
        Energy grid from input file (MeV)
    xs_data : ndarray
        Cross-section data (millibarns)
    reduced_mass : float
        Reduced mass (MeV/c^2)
    
    Returns
    -------
    tuple
        (maxwellian_avg, reaction_rate)
        maxwellian_avg in millibarns
        reaction_rate in cm^3 s^-1 mol^-1
    """
    incident_energy_mev = incident_energy_kev / 1000.0  # Convert keV to MeV
    
    # Neutron velocity (m/s)
    v_t = np.sqrt(2 * incident_energy_mev * C**2 / reduced_mass)
    
    # Maxwellian normalization constant
    const = 2 / (np.sqrt(PI) * incident_energy_mev**2)
    
    # Maxwellian-weighted cross-section integrand
    # σ(E) * E * exp(-E/kT) where kT = incident_energy_mev
    integrand = const * xs_data * energy_data * np.exp(-energy_data / incident_energy_mev)
    
    # Simpson's rule integration
    maxwellian_avg = simpson(integrand, energy_data)
    
    # Reaction rate: <σv> = mxw * v_t * N_A * conversion_factor
    # 1e-25 accounts for: mb->cm² (10^-27) * m/s->cm/s (10^2)
    reaction_rate = maxwellian_avg * v_t * AVO * 1e-25
    
    return maxwellian_avg, reaction_rate


def run_calculations(target_mass, energy_data, xs_data, start_e=5, end_e=100, step_e=5):
    """
    Run reaction rate calculations over energy range.
    
    Parameters
    ----------
    target_mass : float
        Target nucleus mass number
    energy_data : ndarray
        Energy grid (MeV)
    xs_data : ndarray
        Cross-section data (mb)
    start_e : int
        Starting energy (keV)
    end_e : int
        Ending energy (keV)
    step_e : int
        Energy step (keV)
    
    Returns
    -------
    dict
        Dictionary with arrays of results
    """
    reduced_mass = calculate_reduced_mass(target_mass)
    
    energies_kev = np.arange(start_e, end_e + step_e, step_e, dtype=float)
    temperatures_gk = 1.1605 * energies_kev / 100
    maxwellian_avgs = np.zeros_like(energies_kev)
    reaction_rates = np.zeros_like(energies_kev)
    
    for idx, e_kev in enumerate(energies_kev):
        mxw, rr = calculate_reaction_rate(e_kev, energy_data, xs_data, reduced_mass)
        maxwellian_avgs[idx] = mxw
        reaction_rates[idx] = rr
    
    return {
        'energy_kev': energies_kev,
        'temperature_gk': temperatures_gk,
        'maxwellian_mb': maxwellian_avgs,
        'reaction_rate': reaction_rates
    }


def save_results(results, filename):
    """
    Save results to text file with formatted table.
    
    Parameters
    ----------
    results : dict
        Dictionary with result arrays
    filename : str
        Output filename
    """
    with open(filename, 'w') as f:
        f.write("ICARUS v1.22 - Nuclear Reaction Rate Results\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"{'Energy (keV)':>12} {'Temp (GK)':>12} {'<σv> (mb)':>15} {'Rate (cm³/s/mol)':>18}\n")
        f.write("-" * 70 + "\n")
        for i in range(len(results['energy_kev'])):
            f.write(f"{results['energy_kev'][i]:>12.1f} "
                    f"{results['temperature_gk'][i]:>12.4f} "
                    f"{results['maxwellian_mb'][i]:>15.6e} "
                    f"{results['reaction_rate'][i]:>18.6e}\n")
        f.write("-" * 70 + "\n")
        f.write("\nICARUS congratulates you on performing a successful calculation.\n")
    print(f"\nResults saved to: {filename}")


def create_plots(energy_data, xs_data, results):
    """
    Create interactive visualization plots.
    
    Parameters
    ----------
    energy_data : ndarray
        Input energy grid (MeV)
    xs_data : ndarray
        Input cross-section (mb)
    results : dict
        Calculation results
    """
    plt.style.use('default')
    
    # Figure with 2 subplots (side by side)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('ICARUS v1.22 - Nuclear Reaction Analysis', fontsize=14, fontweight='bold')
    
    # Plot 1: Maxwellian-averaged cross-section
    ax1 = axes[0]
    ax1.plot(results['energy_kev'], results['maxwellian_mb'], 'g-s', linewidth=2, markersize=6)
    ax1.set_xlabel('Energy (keV)')
    ax1.set_ylabel('Maxwellian-Averaged σ (mb)')
    ax1.set_title('Maxwellian-Averaged Cross-Section')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Reaction rate vs energy
    ax2 = axes[1]
    ax2.plot(results['energy_kev'], results['reaction_rate'], 'r-o', linewidth=2, markersize=6)
    ax2.set_xlabel('Energy (keV)')
    ax2.set_ylabel('Reaction Rate (cm³ s⁻¹ mol⁻¹)')
    ax2.set_title('Reaction Rate vs Incident Energy')
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='ICARUS v1.22 - Nuclear Reaction Rate Calculator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  python icarus.py ncap.txt --target 152
  python icarus.py ncap.txt --target 152 --output results.txt --plot
        '''
    )
    
    parser.add_argument('input_file', type=str, help='Input cross-section data file')
    parser.add_argument('--target', '-t', type=float, required=True, 
                        help='Target nucleus mass number (A)')
    parser.add_argument('--output', '-o', type=str, default='results.txt',
                        help='Output text filename (default: results.txt)')
    parser.add_argument('--plot', '-p', action='store_true',
                        help='Generate interactive plots')
    parser.add_argument('--start-e', type=int, default=5,
                        help='Starting energy in keV (default: 5)')
    parser.add_argument('--end-e', type=int, default=100,
                        help='Ending energy in keV (default: 100)')
    parser.add_argument('--step-e', type=int, default=5,
                        help='Energy step in keV (default: 5)')
    
    args = parser.parse_args()
    
    # Print logo
    print_logo()
    
    # Load data
    try:
        energy_data, xs_data, n_points = load_cross_section_data(args.input_file)
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    
    # Print input summary
    print(f"A = {args.target}")
    print(f"Input file: {args.input_file}")
    print(f"Data points: {n_points}")
    print("")
    
    # Run calculations
    results = run_calculations(
        target_mass=args.target,
        energy_data=energy_data,
        xs_data=xs_data,
        start_e=args.start_e,
        end_e=args.end_e,
        step_e=args.step_e
    )
    
    # Print results table
    print(f"{'Energy (keV)':>12} {'Temp (GK)':>12} {'<σv> (mb)':>15} {'Rate (cm³/s/mol)':>18}")
    print("-" * 60)
    for i in range(len(results['energy_kev'])):
        print(f"{results['energy_kev'][i]:>12.1f} "
              f"{results['temperature_gk'][i]:>12.4f} "
              f"{results['maxwellian_mb'][i]:>15.6e} "
              f"{results['reaction_rate'][i]:>18.6e}")
    
    # Save results
    save_results(results, args.output)
    
    # Generate plots if requested
    if args.plot:
        create_plots(energy_data, xs_data, results)
    
    # Success message
    print("")
    print("ICARUS congratulates you on performing a successful calculation.")


if __name__ == '__main__':
    main()
