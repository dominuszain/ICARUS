# ICARUS

**Maxwellian-averaged cross sections and thermonuclear reaction rates**

[![DOI](https://img.shields.io/badge/DOI--blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)

ICARUS computes Maxwellian-averaged cross sections (MACS) and thermonuclear
reaction rates for neutron-capture reactions. It takes tabulated cross sections
as input and integrates them over a Maxwellian energy distribution. The code was
developed for my bachelor's thesis on heavy nuclei and can be coupled
with reaction-network codes such as YAHFC.

**Author:** Zain Ul Abideen
**Bachelor's thesis:** [ICARUS]() (Institute of Space Technology, Islamabad, 2024)
**GitHub repository:** [dominuszain/ICARUS](https://github.com/dominuszain/ICARUS)

## Latest version

The current implementation is the Python version in [`v1.22/`](v1.22/). The
earlier Fortran and Scilab implementations are kept in [`legacy/`](legacy/) for
reference.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r v1.22/requirements.txt
cd v1.22
python icarus.py ncap.txt --target 152
```

This writes `results.txt` with four columns: incident neutron energy (keV),
equivalent stellar temperature (GK), Maxwellian-averaged cross section (mb), and
reaction rate (cm³ s⁻¹ mol⁻¹).

Useful options:

```bash
python icarus.py ncap.txt --target 152 --plot
python icarus.py ncap.txt --target 87 --output rb87_results.txt
python icarus.py ncap.txt --target 152 --start-e 1 --end-e 200 --step-e 10
```

See [`v1.22/README.md`](v1.22/README.md) for the complete CLI reference, input
format, and output description. [`v1.22/DETAIL.md`](v1.22/DETAIL.md) documents the
implementation in more depth.

## Input format

The input file has two whitespace-separated columns:

| Column | Quantity | Unit |
|---|---|---|
| 1 | Incident neutron energy | MeV |
| 2 | Cross section | millibarn |

Sample inputs are in [`v1.22/ncap.txt`](v1.22/ncap.txt),
[`v1.22/sample_input.txt`](v1.22/sample_input.txt), and
[`v1.22/152gd.txt`](v1.22/152gd.txt).

## Method

The Maxwellian-averaged cross section is evaluated as

```
<σv> = (2 / (√π (kT)²)) ∫ σ(E) E exp(-E/kT) dE
```

The integral uses Simpson's rule, and the upper limit is replaced by a
sufficiently large finite energy. The reaction rate follows from `<σv>`, the
neutron velocity, and Avogadro's number. The v1.22 implementation also detects
the number of input points automatically and can produce plots of the cross
section and reaction rate.

The legacy Fortran and Scilab versions in `legacy/` use a Riemann-sum
approximation with a finite cutoff; v1.22 uses Simpson's rule.

## Repository layout

```
.
├── v1.22/                  # current Python implementation
│   ├── icarus.py
│   ├── README.md
│   ├── DETAIL.md
│   ├── requirements.txt
│   └── sample inputs
├── legacy/                 # earlier Fortran and Scilab versions
│   ├── fortran-v1.21/
│   ├── fortran-prototypes/
│   ├── fortran-small-update/
│   └── scilab/
├── docs/                   # reference plots
├── CITATION.cff
└── LICENSE
```

## Citation

If you use ICARUS, please cite the thesis and the code:

```bibtex
@misc{ulabideen_icarus_thesis_2024,
  author = {Ul Abideen, Zain},
  title  = {ICARUS},
  year   = {2024},
  note   = {Bachelor's thesis, Institute of Space Technology, Islamabad},
  doi    = {},
  url    = {}
}

@software{ulabideen_icarus_2026,
  author = {Ul Abideen, Zain},
  title  = {ICARUS: Maxwellian-averaged cross sections and thermonuclear reaction rates},
  year   = {2026},
  url    = {https://github.com/dominuszain/ICARUS}
}
```

## License

The code is released under the [MIT License](LICENSE).
