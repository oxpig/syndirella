<div align="center">
  <img src="logos/Full.png" alt="Syndirella Logo" width="400"/>
</div>

![GitHub Tag](https://img.shields.io/github/v/tag/xchem/syndirella?include_prereleases&label=PyPI&link=https%3A%2F%2Fpypi.org%2Fproject%2Fxchem-syndirella%2F)
![release](https://img.shields.io/github/actions/workflow/status/xchem/syndirella/release.yaml?label=release&link=https%3A%2F%2Fgithub.com%2Fxchem%2Fsyndirella%2Factions%2Fworkflows%2Frelease.yaml)

Syndirella (Synthesis Directed Elaborations) is a tool for generating and scoring synthetically practical elaborations of molecules designed from fragment screens. Starting from the retrosynthetic analysis of fragment merged molecules, superstructures of the original reactants are found by catalog search and filtered by selectivity issues. The elaborated final products are defined by reaction SMIRKS patterns and energy minimized in the protein with restraints to experimental data.

**Preprint:** [https://chemrxiv.org/engage/chemrxiv/article-details/68d4f08ff416303770403a44](https://chemrxiv.org/engage/chemrxiv/article-details/68d4f08ff416303770403a44)

## Documentation

**[Full Documentation](https://syndirella.readthedocs.io/)** - Complete user guide, API reference, and examples

[![Documentation Status](https://readthedocs.org/projects/syndirella/badge/?version=latest)](https://syndirella.readthedocs.io/en/latest/?badge=latest)

## Quick Start

### Installation

```bash
conda create -n syndirella python=3.10
conda activate syndirella
pip install syndirella
pip install aizynthfinder
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'  # run once, for Fragmenstein
```

**Note:** PyRosetta is required by Fragmenstein to operate. Syndirella does not use PyRosetta by default to score compounds. PyRosetta is available for academic and non-commercial use (see [PyRosetta License](https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.PyRosetta.md)).

### Troubleshooting:

**cgrtools fails to install:**
```bash
pip install "cython<3.2"
conda install -c conda-forge c-compiler cxx-compiler
pip install --no-build-isolation cgrtools
```

### Install with [HIPPO](https://github.com/mwinokan/HIPPO) (optional)

```bash
pip install aizynthfinder --no-deps # HIPPO requires `rdkit>=2024.9.6`; AiZynthFinder declares `rdkit<2024`. Install in this order so pip does not pull AiZynthFinder’s RDKit
pip install syndirella[hippo]
conda install -c conda-forge chemicalite
```
See the [HIPPO documentation](https://hippo-docs.winokan.com) for full setup and version compatibility.

### Basic Usage

1. **Setup AiZynthFinder** (required for retrosynthesis):

   ⚠️ WARNING: This step downloads large model files (~750MB total) if run the first time. These are required to run AiZynthFinder.

   Files will be downloaded to: `[syndirella_package_path]/aizynth/`

   ```bash
   syndirella setup-aizynth
   ```

2. **Run the pipeline**:
   ```bash
   syndirella run --input /absolute/path/to/input.csv --output /absolute/path/to/output_directory --templates /absolute/path/to/templates/ --hits_path /absolute/path/to/fragments.sdf
   ```

   **Note**: All paths must be absolute paths.

### Example Input

Create a CSV file with your fragment data. Hit names must match the exact names in your SDF file:

```csv
smiles,hit1,hit2,hit3,template,compound_set
O=C(NC1CC(C(F)(F)F)C1)c1cc2ccsc2[nH]1,A71EV2A-x0556_A_147_1_A71EV2A-x0526+A+147+1,A71EV2A-x0566_A_147_1_A71EV2A-x0526+A+147+1,,Ax0310a_apo-desolv,my_compound_set
CC(=O)Nc1cc(CC(=O)NCC(NC(=O)CCl)c2cccnc2)c(NC(C)=O)nn1,A71EV2A-x1346_A_250_1_A71EV2A-x0526+A+147+1,,,Ax1346a_apo-desolv,my_compound_set
```

**For complete examples and templates, see the [run_syndirella_example](examples/run_syndirella_example/) directory.**

### Output Files

The main output file to read is the **`structured_output.pkl.gz`** file, which contains:
- Complete synthesis routes with reaction details
- Placement information with energy values
- Structural quality metrics
- Product information and success flags

This file contains all the information needed to analyze and reproduce the elaborations.

## Requirements

- Python 3.10
- RDKit
- AiZynthFinder
- Fragmenstein
- Additional dependencies (see [pyproject.toml](pyproject.toml))

## Learn More

- **[Examples](examples/)** - Jupyter notebooks and sample data
- **[User Guide](https://syndirella.readthedocs.io/en/latest/user_guide_w_aizynth.html)** - Detailed usage instructions
- **[Reaction Constants](https://syndirella.readthedocs.io/en/latest/configuration.html)** - Advanced reaction settings and options

## License

This project is licensed under the MIT License.
