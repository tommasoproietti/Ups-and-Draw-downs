# Reproducibility Package for: Ups and (Draw)downs

**Date assembled:** 2025-07-04  
**Author: Tommaso Proietti – University of Rome Tor Vergata
          Email: tommaso.proietti@uniroma2.it 
 
 ---

## 📁 Repository Structure

- `data/`: MATLAB `.mat` files and CSVs used for analysis and forecasting.
- `scripts/`: MATLAB scripts for generating all figures and tables in the paper.
- `functions/`: Custom MATLAB functions used in the scripts.
- `README.md`: This document.
- `LICENSE`: License information (e.g., MIT, CC-BY-NC).
- `reproduce_all.m`: Master script to reproduce all figures and tables (optional).

---

## 💻 Computing Environment

- **Language**: MATLAB R2023a  
- **Toolboxes used**: Statistics and Machine Learning, Econometrics Toolbox, Parallel Computing  
- **Platform**: Windows 
- **Expected runtime**: ~ 20 minutes on Intel i7, 16 GB RAM
- **Special setup**:  parallel computing enhances speed

---

## 🔧 Running the Reproducibility Check

Each script is self-contained and will generate the corresponding tables/figures as in the paper. The typical workflow is:

### To reproduce Section 2:

```matlab
run('scripts/sSP500_Figure1Table1.m')
