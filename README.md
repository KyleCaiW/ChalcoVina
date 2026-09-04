ChalcoVina: Enabling Docking for Chalcogen Bond
----------------------------------------------------

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://img.shields.io/static/v1?label=Paper&message=10.1021/acs.jcim.6c02680&color=blue)](https://doi.org/10.1021/acs.jcim.6c02680)

**ChalcoVina** is an augmented version of the AutoDock Vina 1.2.7, specifically developed to simultaneously evaluate intermolecular and intramolecular chalcogen bonds. For more details, please refer to the paper.

## Installation

You can either download the **Linux executable** from the Releases page, or compile it from source (same as AutoDock Vina; see: [readthedocs.org](https://autodock-vina.readthedocs.io/en/latest/)).

## Quick Start

ChalcoVina maintains command-line compatibility with Vina. You can run it using:

```
ChalcoVina --receptor protein.pdbqt --ligand ligand.pdbqt --config config.txt --out output.pdbqt --chbmap ChBGrid.map
```

## Citations

If you use ChalcoVina in your research, please cite the following papers:

* [W. Cai, Z. Li, Y. Wang, W. Yu, Y. Sun, Y. Xu, Q. Zhu, and Y. Zou. (2026). ChalcoVina: Enhancing Molecular Docking Accuracy by Incorporating Chalcogen Bond-Aware Scoring. Journal of Chemical Information and Modeling.](https://pubs.acs.org/doi/10.1021/acs.jcim.6c02680)
* [J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00203)
