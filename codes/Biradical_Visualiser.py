import os

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
from stmol import showmol, makeobj, add_hover


def plot_2d_struct_span_ketcher(smiles_2dfn):
    mol = Chem.MolFromSmiles( smiles_2dfn )
    mol = Chem.AddHs ( mol )
    if mol:
        rdDepictor.Compute2DCoords ( mol )  # Correct usage
        mol = Chem.RemoveHs ( mol )
        img = Draw.MolToImage(mol,size=(400,400))
        st.image(img, caption=f"{biradical_chosen}", use_container_width=True)
        st.subheader("You can edit the structure of the molecule here")
        st.write('Credit: [Streamlit-Ketcher](https://github.com/mik-laj/streamlit-ketcher)')
        st_ketcher( smiles_2dfn )


def generate_3d_molecule(smiles_3dfn):
    try:
        mol = Chem.MolFromSmiles( smiles_3dfn )
        mol3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3 ()
        params.randomSeed = 0xf00d  # optional random seed for reproducibility
        AllChem.EmbedMolecule ( mol3d , params )
        AllChem.MMFFOptimizeMolecule ( mol3d )
        # Convert to MolBlock
        mol_block = Chem.MolToMolBlock(mol3d)
        py3dmol_obj_retn = makeobj ( mol_block , molformat='mol' , style='stick' , background='white' )
        add_hover ( py3dmol_obj_retn , backgroundColor='white' , fontColor='black' )
        return py3dmol_obj_retn
    except Exception as e:
        st.error(f"Error generating 3D molecule: {e}")
        return None

# Load SMILES data

smile_file = os.path.normpath(os.getcwd() + os.sep + os.pardir + '/dep/smilecodes_biradicals.xlsx')
df_smile_file = pd.read_excel(smile_file)
df_smile_file = df_smile_file.sort_values(by="Biradical Name")
biradical_choice = df_smile_file["Biradical Name"]

# Streamlit app
st.title("Biradical Structure Viewer")

# Select biradical
biradical_chosen = st.selectbox("Choose a biradical", options=biradical_choice, index=2)
smiles = df_smile_file[df_smile_file["Biradical Name"] == biradical_chosen]['SMILES Code'].values[0]

st.divider()
st.header("2D Structure Viewer")
plot_2d_struct_span_ketcher(smiles)


# 3D structure viewer
st.divider()
st.header("3D Structure Viewer")
st.markdown('''
This is a 3D structure generated from the .mol file by RDKit.
It is not a DFT optimised structure. Please consider it as a representative structure
The 3D view is generated using the package stmol:
Nápoles-Duarte JM, Biswas A, Parker MI, Palomares-Baez JP, Chávez-Rojo MA and Rodríguez-Valdez LM (2022) Stmol: A component for building interactive molecular visualizations within streamlit web-applications. Front. Mol. Biosci. 9:990846. doi: 10.3389/fmolb.2022.990846
''')

py3dmol_obj = generate_3d_molecule(smiles)
showmol(py3dmol_obj, height=500, width=800)




