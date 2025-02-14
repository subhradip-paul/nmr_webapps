import streamlit as st
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from streamlit_ketcher import st_ketcher

# Load SMILES data
smile_file = os.path.normpath(os.getcwd() + os.sep + os.pardir + '/dep/smilecodes_biradicals.xlsx')
df_smile_file = pd.read_excel(smile_file)
df_smile_file = df_smile_file.sort_values(by="Biradical Name")
biradical_choice = df_smile_file["Biradical Name"]

# Streamlit app
st.title("SMILES to 2D Molecule Viewer")

# Select biradical
biradical_chosen = st.selectbox("Choose a biradical", options=biradical_choice, index=2)
smiles = df_smile_file[df_smile_file["Biradical Name"] == biradical_chosen]['SMILES Code'].values[0]
# st.write(f"SMILES Code: {smiles}")

if smiles:
    try:
        # Generate molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)

        # mol = Chem.SanitizeMol(mol)
        if mol:
            # Draw molecule
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img, caption=f"{biradical_chosen}")
            st_ketcher(smiles)
        else:
            st.error("Invalid SMILES string. Please enter a valid one.")
    except Exception as e:
        st.error(f"Error generating molecule: {e}")