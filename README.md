# PROTEIN-STRUCTURE-ANALYSIS-USING-1A80.PDB-FILE-
SOURCE FILE : 1A80 PDB FILE
# 🧬 Protein Structure Analysis with BioPython (PDB ID: 1A80)

## 📝 Overview

This Jupyter Notebook performs an advanced structural analysis of a protein, specifically the one corresponding to **PDB ID 1A80** (though the code is generic enough for any PDB file).

The analysis leverages the **Bio.PDB** module within the BioPython library to parse the structure, extract atomic coordinates, calculate physical properties (like residue-residue distances and dihedral angles), and visualize key structural features.

## 🚀 Analysis Workflow

The notebook executes the following key steps using the `Bio.PDB` module, `numpy`, and `matplotlib/seaborn`:

### 1. Setup and Data Loading

* **Libraries**: Installs and imports necessary libraries, including `biopython`, `PDBParser`, `PPBuilder`, `NeighborSearch`, `calc_dihedral`, `numpy`, `matplotlib`, and `seaborn`.
* **File Upload**: Includes code using `google.colab.files.upload()` to allow the user to easily upload the required PDB file (e.g., `1a80.pdb`).
* **Structure Parsing**: A `PDBParser` is initialized to read the uploaded PDB file and create a Bio.PDB `Structure` object.

### 2. Structural Inspection and Calculation

* **Hierarchical Traversal**: The script iterates through the **Model**, **Chain** (Chain ID: 'A'), and **Residue** levels of the structure, printing the name and index of each residue (including non-amino acid residues like **NDP** and **HOH**).
* **Cα Atom Extraction**: All **Cα (Alpha Carbon)** atoms belonging to standard amino acids are extracted into a list (`ca_atoms`).
* **Residue Neighbors**: A `NeighborSearch` object is used to find all atoms within a **5.0 Å** radius of the first Cα atom (Residue 2, THR), demonstrating basic neighbor analysis.

### 3. Visualization of Backbone Geometry

* **Cα-Cα Distance Plot**:
    * The **distance between consecutive Cα atoms** is calculated and plotted against the residue index.
    * This plot is useful for quickly identifying potential gaps, missing residues, or large conformational changes in the protein backbone, as the distance between adjacent Cα atoms in a peptide bond is typically very consistent (around 3.8 Å).
* **Ramachandran Plot (φ vs ψ)**:
    * The **Phi ($\phi$) and Psi ($\psi$) backbone dihedral angles** are calculated for every amino acid residue (excluding termini).
    * The resulting angle pairs are plotted in a **Ramachandran plot**, showing the allowed and favored regions of protein secondary structure (alpha-helices and beta-sheets).

### 4. Sequence and Hydrophobicity Mapping

* **Amino Acid Frequency**: A **bar chart** is generated showing the count of each type of amino acid residue present in the protein structure.
* **Hydrophobic Residue Mapping**:
    * A set of standard **hydrophobic residues** (`ALA, VAL, LEU, MET, PHE, TRP, PRO`) is defined.
    * The Cα coordinates of these hydrophobic residues are extracted.
    * A **2D scatter plot** is generated to map the spatial locations of these residues in the X-Y plane, helping to visualize the hydrophobic core or patches on the protein surface.

### 5. Contact Map / Distance Heatmap

* **Cα-Cα Distance Matrix**: A complete distance matrix is computed between all Cα atoms in the structure using `scipy.spatial.distance.pdist`.
* **Heatmap Visualization**: A **heatmap** is generated using `seaborn` to display the distance matrix. This is effectively a **contact map**, where dark regions (low values) indicate residues that are close in 3D space, regardless of their position in the primary sequence.

## 💾 Data Requirement

The notebook requires the user to upload a single PDB format file:

* **Required File**: A valid PDB file (e.g., `1a80.pdb`).

## 🛠️ Key Libraries

* `Bio.PDB`: For parsing, manipulating, and analyzing protein structure files.
* `numpy`: For numerical operations, especially with coordinates and angle calculations.
* `matplotlib` & `seaborn`: For visualizing the backbone trace, Cα distances, residue frequency, and the distance heatmap.
* `scipy.spatial.distance`: Used for generating the distance matrix efficiently.
