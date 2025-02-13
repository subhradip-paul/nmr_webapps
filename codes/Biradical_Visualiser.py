import streamlit as st
from stmol import showmol
import py3Dmol
from rdkit import Chem
import os

# Streamlit App
st.title("Molecule Viewer from MOL File")

# Specify the file path
script_dir = os.path.dirname(__file__)
biradical_choice = ["btbk", "bturea", "bctbk", "abk", "abu", "bmtbk", "tekpol",
                    "pypol", "amupol", "tekpolcbm", "pypolcbm","amupolcbm",
                    "asympolpok", "asympoltek","hydropol", "tekpolcho", "naphpol",
                    "cbctol"]
biradical_chosen = st.selectbox("Choose a biradical", options=biradical_choice, index=2)
st.write(script_dir, os.pardir)
mol_file_path = os.path.join(os.pardir, f'/molfiles_biradicals/{biradical_chosen}.mol')

try:
    # Read MOL file content
    with open(mol_file_path, "r") as f:
        mol_data = f.read()

    # Convert to RDKit Molecule
    mol = Chem.MolFromMolBlock(mol_data)

    if mol:
        # Convert to 3Dmol.js format
        mol_block = Chem.MolToMolBlock(mol)

        # Create 3Dmol.js view
        viewer = py3Dmol.view(width=500, height=500)
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.3}})
        viewer.zoomTo()

        # Show the molecule in Streamlit using stmol
        showmol(viewer, height=500, width=500)

    else:
        st.error("Invalid MOL file. Please check the file.")

except FileNotFoundError:
    st.error(f"File not found: {mol_file_path}")
except Exception as e:
    st.error(f"An error occurred: {e}")