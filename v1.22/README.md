# ICARUS v1.22

**Nuclear Reaction Rate Calculator**

A Python tool for calculating Maxwellian-averaged cross-sections and nuclear reaction rates for neutron capture reactions.

---

## Overview

ICARUS calculates thermonuclear reaction rates for neutron-induced reactions using input cross-section data. It integrates the cross-section over a Maxwellian energy distribution to produce temperature-dependent reaction rates.

**Original Author:** Zain Ul Abideen (FYP, BS-SS09-IST)  
**Contact:** dominuszain@gmail.com  
**Version:** 1.22 (Python implementation)  
**Modernized with:** Opencode + Qwen3.5 397B

---

## Features

- Maxwellian-averaged cross-section calculation using Simpson's rule integration
- Reaction rate computation in cm³ s⁻¹ mol⁻¹
- Configurable energy range and temperature grid
- Formatted text output with publication-ready tables
- Interactive visualization plots (optional)
- Automatic detection of input data points

---

## Requirements

- Python 3.8+
- NumPy
- SciPy
- Matplotlib

---

## Installation

### 1. Clone or download the repository

```bash
cd /path/to/icarus
```

### 2. Create a virtual environment (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate  # On Linux/macOS
# or
.venv\Scripts\activate     # On Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install numpy scipy matplotlib
```

---

## Usage

### Basic Usage

```bash
python icarus.py <input_file> --target <A>
```

**Example:**

```bash
python icarus.py ncap.txt --target 152
```

### Command-Line Arguments

| Argument | Short | Required | Default | Description |
|----------|-------|----------|---------|-------------|
| `input_file` | – | Yes | – | Input cross-section data file |
| `--target` | `-t` | Yes | – | Target nucleus mass number (A) |
| `--output` | `-o` | No | `results.txt` | Output file name |
| `--plot` | `-p` | No | False | Generate interactive plots |
| `--start-e` | – | No | 5 | Starting energy in keV |
| `--end-e` | – | No | 100 | Ending energy in keV |
| `--step-e` | – | No | 5 | Energy step in keV |

### Examples

**1. Standard calculation:**

```bash
python icarus.py ncap.txt --target 152
```

**2. With plots:**

```bash
python icarus.py ncap.txt --target 152 --plot
```

**3. Custom output file:**

```bash
python icarus.py ncap.txt --target 87 --output rb87_results.txt
```

**4. Extended energy range:**

```bash
python icarus.py ncap.txt --target 152 --start-e 1 --end-e 200 --step-e 10
```

---

## Input File Format

The input file should contain two columns:

| Column 1 | Column 2 |
|----------|----------|
| Energy (MeV) | Cross-section (millibarns) |

**Example (`ncap.txt`):**

```
0.001       4539.385
0.002       2422.342
0.003       1692.983
0.004       1322.956
0.005       1099.609
...
```

**Notes:**
- Energy values must be in **MeV**
- Cross-section values must be in **millibarns (mb)**
- Data points should cover the energy range of interest
- File can have any number of data points (auto-detected)

---

## Output

### Text Output (`results.txt`)

The default output is a formatted text file containing:

```
ICARUS v1.22 - Nuclear Reaction Rate Results
======================================================================

Energy (keV)    Temp (GK)       <σv> (mb)   Rate (cm³/s/mol)
----------------------------------------------------------------------
         5.0       0.0580    1.051638e+03       6.212389e+07
        10.0       0.1161    6.911345e+02       5.773908e+07
        ...
----------------------------------------------------------------------

ICARUS congratulates you on performing a successful calculation.
```

| Column | Quantity | Units |
|--------|----------|-------|
| Energy | Incident neutron energy | keV |
| Temp | Equivalent temperature | GK (gigakelvin) |
| <σv> | Maxwellian-averaged cross-section | millibarns |
| Rate | Reaction rate | cm³ s⁻¹ mol⁻¹ |

### Plots (with `--plot` flag)

Two interactive plots are generated:

1. **Maxwellian-Averaged Cross-Section** vs Energy
2. **Reaction Rate** vs Energy (logarithmic y-axis)

---

## Physics Background

### Reaction Rate Formula

The thermonuclear reaction rate is calculated as:

```
rr = <σv> × v × N_A × 10⁻²⁵
```

Where:
- `<σv>` = Maxwellian-averaged cross-section (mb)
- `v` = neutron velocity (m/s)
- `N_A` = Avogadro's number (6.022 × 10²³ mol⁻¹)
- `10⁻²⁵` = unit conversion factor (mb→cm² × m/s→cm/s)

### Maxwellian Average

The Maxwellian-averaged cross-section is computed by integrating:

```
<σv> = ∫ σ(E) × E × exp(-E/kT) dE × [2 / (√π × (kT)²)]
```

Using Simpson's rule for numerical integration.

### Temperature Conversion

```
T (GK) = 1.1605 × E (keV) / 100
```

---

## Version History

| Version | Changes |
|---------|---------|
| 1.21 | Original Fortran code |
| 1.22 | Python implementation with Simpson's rule integration, auto-detection of data points, text output, improved plots |

---

## License

This code was developed as part of a Final Year Project (FYP). Contact the author for licensing information.

---

## Troubleshooting

### Common Issues

**1. "File not found" error:**
- Ensure the input file path is correct
- Use absolute path if needed

**2. Import errors:**
- Activate the virtual environment: `source .venv/bin/activate`
- Install dependencies: `pip install numpy scipy matplotlib`

**3. Plot window doesn't appear:**
- Ensure you have a display server running
- Try saving plots instead (modification required)

---

## Acknowledgments

Original Fortran implementation by Zain Ul Abideen.  
Python upgrade maintains full compatibility while improving accuracy and usability.
