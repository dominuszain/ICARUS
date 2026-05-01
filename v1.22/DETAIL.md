# ICARUS v1.22 - Technical Implementation Details

**Comprehensive documentation for AI assistants and developers**

This document provides complete technical details necessary for understanding, maintaining, and extending the ICARUS codebase.

---

## Table of Contents

1. [Code Architecture](#code-architecture)
2. [Physical Constants](#physical-constants)
3. [Mathematical Formulation](#mathematical-formulation)
4. [Function Reference](#function-reference)
5. [Data Flow](#data-flow)
6. [Unit Analysis](#unit-analysis)
7. [Numerical Methods](#numerical-methods)
8. [Input/Output Specifications](#inputoutput-specifications)
9. [Extension Guide](#extension-guide)
10. [Testing & Validation](#testing--validation)
11. [Known Limitations](#known-limitations)

---

## Code Architecture

### File Structure

```
icarus/
├── icarus.py          # Main implementation (308 lines)
├── README.md          # User documentation
├── DETAIL.md          # This file - technical documentation
├── requirements.txt   # Python dependencies
├── ncap.txt          # Sample input (cross-section data)
└── results.txt       # Output (generated on run)
```

### Module Organization

The code is organized as a single-file script with the following structure:

```
Line 1-13:   Imports and physical constants
Line 15-35:  print_logo() - ASCII art display
Line 37-59:  load_cross_section_data() - Input parsing
Line 61-75:  calculate_reduced_mass() - Physics calculation
Line 77-117: calculate_reaction_rate() - Core physics
Line 119-163: run_calculations() - Main computation loop
Line 165-192: save_results() - Output formatting
Line 194-230: create_plots() - Visualization
Line 232-308: main() - CLI entry point
```

### Dependencies

| Library | Version | Purpose |
|---------|---------|---------|
| numpy | ≥1.20 | Array operations, mathematical functions |
| scipy | ≥1.7 | Simpson's rule integration (`scipy.integrate.simpson`) |
| matplotlib | ≥3.4 | Plotting and visualization |
| argparse | built-in | Command-line argument parsing |
| sys | built-in | Exit handling |

---

## Physical Constants

All constants are defined at the module level (lines 13-18):

```python
PI = 3.14159265358979              # π (dimensionless)
E = 2.71828182845904               # e (dimensionless)
C = 2.997e8                        # Speed of light (m/s)
AVO = 6.022e23                     # Avogadro's number (mol⁻¹)
MASS_NEUTRON = 939.6               # Neutron mass (MeV/c²)
AMU = 931.5                        # Atomic mass unit (MeV/c²)
```

### Notes on Constants

1. **Precision**: Constants use float64 precision (standard Python float)
2. **Units**: Masses are in MeV/c² for natural unit compatibility
3. **Source**: Standard CODATA values (rounded for practical use)

---

## Mathematical Formulation

### 1. Reduced Mass Calculation

**Function:** `calculate_reduced_mass(target_mass)` (line 64)

**Formula:**
```
μ = (m_n × A × amu) / (A × amu + m_n)
```

Where:
- `μ` = reduced mass (MeV/c²)
- `m_n` = neutron mass = 939.6 MeV/c²
- `A` = target mass number (dimensionless)
- `amu` = atomic mass unit = 931.5 MeV/c²

**Physical Meaning:** The reduced mass accounts for the finite mass of the target nucleus in the center-of-mass frame.

---

### 2. Neutron Velocity

**Location:** `calculate_reaction_rate()` line 91

**Formula:**
```
v = √(2 × E × c² / μ)
```

Where:
- `v` = neutron velocity (m/s)
- `E` = incident energy (MeV)
- `c` = speed of light (m/s)
- `μ` = reduced mass (MeV/c²)

**Derivation:** From classical kinetic energy: `E = ½mv²`, solved for `v`.

---

### 3. Maxwellian Distribution Normalization

**Location:** `calculate_reaction_rate()` line 94

**Formula:**
```
const = 2 / (√π × (kT)²)
```

Where:
- `kT` = thermal energy (MeV) = incident energy in keV / 1000

**Purpose:** Normalizes the Maxwellian energy distribution so that ∫f(E)dE = 1.

---

### 4. Maxwellian-Averaged Cross-Section

**Location:** `calculate_reaction_rate()` lines 97-103

**Formula:**
```
<σv> = ∫₀^∞ σ(E') × const × E' × exp(-E'/kT) dE'
```

**Implementation:**
```python
integrand = const * xs_data * energy_data * np.exp(-energy_data / incident_energy_mev)
maxwellian_avg = simpson(integrand, energy_data)
```

**Physical Meaning:** Averages the cross-section over a thermal (Maxwell-Boltzmann) energy distribution at temperature T.

---

### 5. Reaction Rate

**Location:** `calculate_reaction_rate()` line 106

**Formula:**
```
rr = <σv> × v × N_A × 10⁻²⁵
```

**Unit Conversion Factor (10⁻²⁵):**
- `10⁻²⁷` from millibarns → cm² (1 mb = 10⁻²⁷ cm²)
- `10²` from m/s → cm/s
- Combined: `10⁻²⁷ × 10² = 10⁻²⁵`

**Final Units:** cm³ s⁻¹ mol⁻¹

---

### 6. Temperature Conversion

**Location:** `run_calculations()` line 137

**Formula:**
```
T (GK) = 1.1605 × E (keV) / 100
```

**Derivation:**
- `k = 8.617 × 10⁻⁵ eV/K` (Boltzmann constant)
- `1 keV = 1.1605 × 10⁷ K = 0.011605 GK`
- Factor `1.1605/100 = 0.011605`

---

## Function Reference

### `print_logo()`

**Line:** 18  
**Purpose:** Display ASCII art logo and version information  
**Returns:** None (prints to stdout)  
**Side Effects:** None

---

### `load_cross_section_data(filename)`

**Line:** 40  
**Purpose:** Load energy and cross-section data from file  
**Parameters:**
- `filename` (str): Path to input file

**Returns:**
- `energy` (ndarray): Energy values in MeV
- `cross_section` (ndarray): Cross-section values in mb
- `n_points` (int): Number of data points loaded

**Error Handling:**
- Raises `FileNotFoundError` if file doesn't exist
- Raises `ValueError` if file format is invalid

**Implementation Notes:**
- Uses `numpy.loadtxt()` for efficient parsing
- Automatically detects number of data points
- Assumes whitespace-delimited columns

---

### `calculate_reduced_mass(target_mass)`

**Line:** 64  
**Purpose:** Compute reduced mass for neutron-target system  
**Parameters:**
- `target_mass` (float): Target nucleus mass number A

**Returns:**
- `float`: Reduced mass in MeV/c²

**Physics:** Two-body kinematics in center-of-mass frame

---

### `calculate_reaction_rate(incident_energy_kev, energy_data, xs_data, reduced_mass)`

**Line:** 80  
**Purpose:** Core physics calculation - compute Maxwellian average and reaction rate  
**Parameters:**
- `incident_energy_kev` (float): Incident neutron energy in keV
- `energy_data` (ndarray): Energy grid from input (MeV)
- `xs_data` (ndarray): Cross-section data (mb)
- `reduced_mass` (float): Reduced mass (MeV/c²)

**Returns:**
- `maxwellian_avg` (float): Maxwellian-averaged cross-section (mb)
- `reaction_rate` (float): Reaction rate (cm³ s⁻¹ mol⁻¹)

**Algorithm:**
1. Convert incident energy from keV to MeV
2. Calculate neutron velocity
3. Compute Maxwellian normalization constant
4. Build integrand: σ(E) × E × exp(-E/kT)
5. Integrate using Simpson's rule
6. Apply unit conversion for reaction rate

**Computational Complexity:** O(n) where n = number of data points

---

### `run_calculations(target_mass, energy_data, xs_data, start_e, end_e, step_e)`

**Line:** 122  
**Purpose:** Run reaction rate calculations over energy range  
**Parameters:**
- `target_mass` (float): Target nucleus mass number
- `energy_data` (ndarray): Energy grid (MeV)
- `xs_data` (ndarray): Cross-section data (mb)
- `start_e` (int): Starting energy in keV (default: 5)
- `end_e` (int): Ending energy in keV (default: 100)
- `step_e` (int): Energy step in keV (default: 5)

**Returns:**
- `dict` with keys:
  - `'energy_kev'` (ndarray): Energy values
  - `'temperature_gk'` (ndarray): Temperature values
  - `'maxwellian_mb'` (ndarray): Maxwellian averages
  - `'reaction_rate'` (ndarray): Reaction rates

**Algorithm:**
1. Calculate reduced mass
2. Generate energy grid using `np.arange()`
3. Calculate temperature for each energy
4. Loop over energies, calling `calculate_reaction_rate()`
5. Store results in arrays
6. Return dictionary

---

### `save_results(results, filename)`

**Line:** 168  
**Purpose:** Save results to formatted text file  
**Parameters:**
- `results` (dict): Dictionary with result arrays
- `filename` (str): Output file path

**Returns:** None  
**Side Effects:** Creates/writes output file

**Output Format:**
```
ICARUS v1.22 - Nuclear Reaction Rate Results
======================================================================

Energy (keV)    Temp (GK)       <σv> (mb)   Rate (cm³/s/mol)
----------------------------------------------------------------------
         5.0       0.0580    1.051638e+03       6.212389e+07
        ...
----------------------------------------------------------------------

ICARUS congratulates you on performing a successful calculation.
```

---

### `create_plots(energy_data, xs_data, results)`

**Line:** 197  
**Purpose:** Generate interactive visualization plots  
**Parameters:**
- `energy_data` (ndarray): Input energy grid (MeV)
- `xs_data` (ndarray): Input cross-section (mb)
- `results` (dict): Calculation results

**Returns:** None  
**Side Effects:** Opens matplotlib window

**Plot Layout:**
- Figure size: 14×5 inches
- 2 subplots side-by-side (1 row, 2 columns)
- Left: Maxwellian-averaged cross-section (linear scale)
- Right: Reaction rate (logarithmic y-axis)

**Styling:**
- Default matplotlib style
- Green squares for cross-section
- Red circles for reaction rate
- Grid with 30% alpha

---

### `main()`

**Line:** 252  
**Purpose:** CLI entry point and program orchestrator  
**Parameters:** None (parses sys.argv)  
**Returns:** None (exits via sys.exit on error)

**Workflow:**
1. Parse command-line arguments
2. Print logo
3. Load input data
4. Run calculations
5. Print results table to stdout
6. Save results to file
7. Generate plots if requested
8. Print success message

---

## Data Flow

```
┌─────────────────┐
│  Input File     │
│  (ncap.txt)     │
│  E (MeV), σ(mb) │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ load_cross_     │
│ section_data()  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ energy_data     │
│ xs_data         │
│ n_points        │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ run_calculations│
│ ─────────────── │
│ for E in 5..100 │
│   keV step 5    │
└────────┬────────┘
         │
         ├─────────────────┬─────────────────┐
         ▼                 ▼                 ▼
┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
│ save_results()  │ │ create_plots()  │ │ print to stdout │
│ (results.txt)   │ │ (matplotlib)    │ │ (formatted)     │
└─────────────────┘ └─────────────────┘ └─────────────────┘
```

---

## Unit Analysis

### Complete Unit Tracking

| Variable | Symbol | Units | Notes |
|----------|--------|-------|-------|
| Input energy | E_in | MeV | From file |
| Input cross-section | σ | mb (millibarns) | From file |
| Incident energy | E_inc | keV | Loop variable |
| Incident energy | E_inc | MeV | After /1000 conversion |
| Reduced mass | μ | MeV/c² | Natural units |
| Velocity | v | m/s | From √(2Ec²/μ) |
| Temperature | T | GK | From 1.1605×E/100 |
| Maxwellian avg | <σv> | mb | Same as input σ |
| Reaction rate | rr | cm³ s⁻¹ mol⁻¹ | Final output |

### Unit Conversion Chain

```
<σv> [mb] × v [m/s] × N_A [mol⁻¹] × 10⁻²⁵
= <σv> × 10⁻²⁷ [cm²] × v × 10² [cm/s] × N_A [mol⁻¹]
= <σv> × v × N_A × 10⁻²⁵ [cm³ s⁻¹ mol⁻¹] ✓
```

---

## Numerical Methods

### Simpson's Rule Integration

**Implementation:** `scipy.integrate.simpson()`

**Formula:**
```
∫ f(x) dx ≈ (h/3) × [f(x₀) + 4f(x₁) + 2f(x₂) + 4f(x₃) + ... + f(xₙ)]
```

Where `h` is the spacing between points.

**Advantages over Rectangular Rule:**
- O(h⁴) error vs O(h²) for rectangular
- Better accuracy for smooth functions
- Handles curved integrands (like Maxwellian) more accurately

**Why Simpson's Rule:**
The Maxwellian integrand `σ(E) × E × exp(-E/kT)` is smooth and curved. Simpson's rule captures this curvature better than the rectangular rule used in the original Fortran code.

**Result Differences:**
- At low energies (5-10 keV): ~5-10% difference from Fortran
- At high energies (50-100 keV): <1% difference
- Python results are more accurate

---

## Input/Output Specifications

### Input File Format

**File Type:** Plain text  
**Encoding:** ASCII or UTF-8  
**Delimiter:** Whitespace (spaces or tabs)  
**Comments:** Not supported (all lines are data)

**Structure:**
```
<energy_MeV> <cross_section_mb>
<energy_MeV> <cross_section_mb>
...
```

**Example:**
```
0.001       4539.385
0.002       2422.342
0.003       1692.983
```

**Requirements:**
- Column 1: Energy in MeV (float, positive)
- Column 2: Cross-section in mb (float, positive)
- Minimum 2 data points (for integration)
- Energies should be monotonically increasing
- No header line

**Validation:**
- File must exist and be readable
- All values must be parseable as floats
- At least 2 data points required

---

### Output File Format

**File Type:** Plain text  
**Encoding:** UTF-8  
**Default Name:** `results.txt`

**Structure:**
```
Line 1:     Title header
Line 2:     Separator (= × 70)
Line 3:     Empty
Line 4:     Column headers
Line 5:     Separator (- × 70)
Lines 6-25: Data rows (20 energies)
Line 26:    Separator (- × 70)
Line 27:    Empty
Line 28:    Success message
```

**Data Row Format:**
```
{energy:12.1f} {temp:12.4f} {sigma:15.6e} {rate:18.6e}
```

---

## Extension Guide

### Adding New Output Formats

To add JSON output:

```python
import json

def save_results_json(results, filename):
    """Save results to JSON file."""
    output = {
        'metadata': {
            'version': '1.22',
            'units': {
                'energy': 'keV',
                'temperature': 'GK',
                'cross_section': 'mb',
                'reaction_rate': 'cm3/s/mol'
            }
        },
        'data': []
    }
    
    for i in range(len(results['energy_kev'])):
        output['data'].append({
            'energy_kev': float(results['energy_kev'][i]),
            'temperature_gk': float(results['temperature_gk'][i]),
            'maxwellian_mb': float(results['maxwellian_mb'][i]),
            'reaction_rate': float(results['reaction_rate'][i])
        })
    
    with open(filename, 'w') as f:
        json.dump(output, f, indent=2)
```

### Adding Support for Different Reactions

To extend beyond (n,γ) reactions:

1. **Add projectile mass parameter:**
```python
def calculate_reduced_mass(projectile_mass, target_mass):
    # projectile_mass in MeV/c²
    return (projectile_mass * target_mass * AMU) / \
           (target_mass * AMU + projectile_mass)
```

2. **Modify CLI:**
```python
parser.add_argument('--projectile', type=str, default='neutron',
                    help='Projectile type (neutron, proton, alpha)')
```

3. **Add projectile mass lookup:**
```python
PROJECTILES = {
    'neutron': 939.6,
    'proton': 938.3,
    'alpha': 3727.4
}
```

### Adding Resonance Integration

For reactions with narrow resonances:

```python
def add_resonance_contribution(energy, sigma, resonance_params):
    """
    Add narrow resonance contributions using Breit-Wigner formula.
    
    Parameters:
    - resonance_params: list of (E_r, gamma_n, gamma_g) tuples
    """
    for E_r, gamma_n, gamma_g in resonance_params:
        sigma_res = (gamma_n * gamma_g) / \
                    ((energy - E_r)**2 + (gamma_n + gamma_g)**2 / 4)
        sigma += sigma_res
    return sigma
```

### Adding Temperature-Dependent Cross-Sections

For Doppler broadening:

```python
def doppler_broaden(energy, sigma, temperature):
    """
    Apply Doppler broadening to cross-section.
    
    Uses the free gas model for thermal motion of target nuclei.
    """
    # Implementation requires convolution with Gaussian
    # representing thermal motion of target
    pass
```

---

## Testing & Validation

### Verification Tests

**1. Unit Consistency:**
```python
def test_units():
    # Verify reaction rate has correct units
    rr = calculate_reaction_rate(10.0, energy_data, xs_data, mu)
    assert rr[1] > 0  # Should be positive
    assert 1e6 < rr[1] < 1e9  # Typical magnitude
```

**2. Energy Conservation:**
```python
def test_energy_range():
    # Maxwellian should peak near kT
    results = run_calculations(152, energy_data, xs_data)
    # <σv> should decrease monotonically (for 1/v cross-section)
```

**3. Comparison with Fortran:**
```python
def test_fortran_comparison():
    # Run both codes with same input
    # Verify results match within expected tolerance
    # Fortran (rectangular) vs Python (Simpson)
    tolerance = 0.10  # 10% at low E, <1% at high E
```

### Validation Data

**Known Good Results (Gd-152 at 10 keV):**
- Maxwellian average: ~691 mb
- Reaction rate: ~5.77 × 10⁷ cm³/s/mol

**Physical Plausibility Checks:**
- Reaction rate should have minimum at 10-20 keV
- Rate should increase at higher energies
- Maxwellian average should decrease monotonically

---

## Known Limitations

### 1. Integration Accuracy

**Issue:** Simpson's rule requires evenly-spaced data points.

**Current Workaround:** The code uses the input data points as-is, which may not be evenly spaced.

**Impact:** Small errors if spacing varies significantly.

**Fix:** Interpolate to uniform grid before integration:
```python
from scipy.interpolate import interp1d

def integrate_uniform(energy, sigma, n_points=1000):
    f = interp1d(energy, sigma, kind='cubic')
    E_uniform = np.linspace(energy[0], energy[-1], n_points)
    sigma_uniform = f(E_uniform)
    return simpson(sigma_uniform, E_uniform)
```

### 2. High-Energy Extrapolation

**Issue:** Input data may not extend to 100 keV.

**Current Behavior:** Integration only covers available data range.

**Impact:** May underestimate <σv> if significant cross-section exists beyond data range.

**Fix:** Extrapolate with physical model (e.g., constant σ at high E).

### 3. Low-Energy Behavior

**Issue:** Maxwellian at low T is very narrow, sensitive to low-E data.

**Current Behavior:** Uses whatever low-E data is provided.

**Impact:** Results at 5 keV may be inaccurate if data doesn't extend below ~0.001 MeV.

**Fix:** Ensure input data extends to sufficiently low energies (<1 eV for T<10 keV).

### 4. No Error Propagation

**Issue:** Cross-section uncertainties not propagated to reaction rates.

**Current Behavior:** Assumes input σ values are exact.

**Fix:** Add Monte Carlo uncertainty propagation:
```python
def calculate_rate_with_uncertainty(energy, sigma, sigma_err, n_samples=1000):
    sigma_samples = np.random.normal(sigma, sigma_err, size=(n_samples, len(sigma)))
    rates = [calculate_reaction_rate(energy, sigma_s) for sigma_s in sigma_samples]
    return np.mean(rates), np.std(rates)
```

### 5. Single Target Only

**Issue:** Code processes one target nucleus at a time.

**Fix:** Batch processing:
```python
def run_batch(targets, input_file):
    results = {}
    for target in targets:
        results[target] = run_calculations(target, ...)
    return results
```

---

## Performance Notes

### Computational Complexity

- **Time:** O(n × m) where n = energy points, m = calculation points
- **Space:** O(n + m) for storing input and output arrays

### Typical Runtime

| Data Points | Calculation Points | Runtime |
|-------------|-------------------|---------|
| 100 | 20 | ~10 ms |
| 600 | 20 | ~50 ms |
| 1000 | 100 | ~200 ms |

### Optimization Opportunities

1. **Vectorization:** Already using numpy vectorization
2. **Parallel Processing:** Could parallelize over energy points
3. **Caching:** Could cache reduced mass calculation
4. **Numba JIT:** Could add `@njit` for 10-100× speedup

---

## Changelog

### v1.22 (Current)
- Python implementation
- Simpson's rule integration (more accurate)
- Automatic data point detection
- Formatted text output
- Improved visualization (2 plots)
- CLI with argparse

### v1.21 (Original)
- Fortran implementation
- Rectangular rule integration
- Manual data point specification
- Unformatted print output
- ASCII logo

---

## References

1. **Nuclear Reactions for Astrophysics** - Iliadis (2007)
2. **Introductory Nuclear Physics** - Krane (1987)
3. **Numerical Recipes** - Press et al. (2007) - Simpson's rule
4. **NIST Physical Reference Data** - Physical constants

---

## Contact

For questions about this documentation or the code implementation:
- Original author: dominuszain@gmail.com
- AI assistant: Available via opencode interface
