# PDB Transmembrane Protein - Membrane Orientation Scanner

Command-line tool to predict the orientation and position of the lipid membrane relative to a transmembrane protein, from a PDB structure. The program centers the protein, generates a scanning grid of directions on a sphere, and identifies the most favorable membrane axis and position based on hydrophobicity, using one of several scanning methods. The result can be visualized directly in PyMOL.

## Features

- Load a transmembrane protein PDB structure and its solvent accessible surface area (SASA) data (`.rsa` file)
- Automatically center the protein at the origin
- Generate a grid of points to scan possible membrane orientations around the protein
- 4 available scanning methods (vectorized and non-vectorized, 2 different approaches)
- Estimate the optimal membrane axis, hydrophobicity score, and membrane center position along that axis
- Visualize the predicted membrane (hydrophobic planes) around the protein in PyMOL
- Benchmark script to compare the performance of the 4 methods

## Requirements

- Python 3.x
- [PyMOL](https://pymol.org/) (for visualization)
- [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/) (to generate `.rsa` SASA files, if you don't already have them)

## Installation

```bash
git clone https://github.com/Eve-Ast/TMDET_BI.git
cd TMDET_BI
pip install -r requirements.txt
```

The project also relies on the standard library (`argparse`, `math`, `os`, `time`) and its internal modules (`Protein`, `Grid`).

## Project structure

```
.
├── src/
│   ├── main.py          # Main script
│   ├── benchmark.py     # Benchmark script for the methods
│   ├── Protein.py        # Protein class
│   └── Grid.py            # Grid class
├── data/                 # Input PDB files
└── result/               # Input RSA (SASA) files
```

> ⚠️ Adjust this layout if it differs from your actual project organization.

## Usage

### Scan a protein

```bash
python src/main.py <pdb_file> [--method <method>]
```

**Arguments:**
- `pdb_file`: name or path of the transmembrane protein PDB file to analyze (the `.pdb` extension can be omitted)
- `--method`: scanning method to use, among:
  - `method1_vec`
  - `method2_vec` *(default)*
  - `method1_novec`
  - `method2_novec`

**Examples:**

```bash
python src/main.py 1PRN.pdb
python src/main.py 1PRN.pdb --method method1_vec
```

**Output:**
- Best detected membrane orientation axis
- Associated hydrophobicity score
- Estimated membrane center position
- A PyMOL window opens with the visualization

### Compare method performance (benchmark)

```bash
python src/benchmark.py <pdb_file> <rsa_file> <n_points>
```

**Example:**

```bash
python src/benchmark.py 1PRN.pdb 1PRN.rsa 1000
```

The script displays, for each of the 4 methods, the execution time and the resulting hydrophobicity score.

## Author

Eve-Angeline STEPHEN

## License

*(to be specified)*
