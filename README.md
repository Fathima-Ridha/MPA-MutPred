# **MPA-MutPred: A Novel Strategy for Accurately Predicting the Binding Affinity Change Upon Mutation in Membrane Protein Complexes**

### **Fathima Ridha**, **M. Michael Gromiha<sup>#</sup>**

---
We developed [MPA-MutPred](https://web.iitm.ac.in/bioinfo2/MPA-MutPred/), a novel method specific for predicting the binding affinity change upon mutation (ΔΔG) in membrane protein-protein complexes.

## Usage Restrictions

The code has been provided for **academic purposes only**.

## Prerequisites

This program is developed using **Python 3** (>v3.8) and requires several Python packages to run. Please ensure the following packages are installed:

- **Bio**
- **pandas**
- **numpy**
- **scipy**
- **scikit-learn**
- **networkx**
- **Generic libraries**: `sys`, `re`, `subprocess`, `os`, `math`,`collections`, ...

### Additional Software Requirements

- **Naccess (v2.1.1)**: Ensure you edit the Naccess executable file with the path of your local installation. You can download it from [Naccess](http://www.bioinf.manchester.ac.uk/naccess/).

- **FoldX (v2024)**: Kindly install the FoldX software from [FoldX](https://foldxsuite.crg.eu/).
- **PSI-BLAST**
- **UniProt db**


## Code Usage

To run the model for a single mutation, execute the following command from within the working directory:

      python3 code_main.py <PDBID> <chain1> <chain2> <mutation>

### Arguments:

- **`<PDBID>`**: The ID of the PDB file (e.g., `1iar`).
- **`<chain1>`**: The identifier for the first chain (e.g., `A`).
- **`<chain2>`**: The identifier for the second chain (e.g., `B`).
- **`<mutation>`**: The mutation string (e.g., `YB13F`).

## Example

Here's an example of how to use the command:

      python3 code_main.py 1iar A B YB13F > output.txt
  
For further queries related to code usage, feel free to reach out via email (bt20d701@smail.iitm.ac.in, gromiha@iitm.ac.in) 


## Citation

Please cite the following article if you use the code in this repository for your research:

**Ridha, F.**, &amp; Gromiha, M. M. (2024). MPA-MutPred: A novel strategy for accurately predicting the binding affinity change upon mutation in membrane protein complexes.
