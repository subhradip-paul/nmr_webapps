import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
from scipy.spatial.transform import Rotation as Rot
from io import StringIO
import py3Dmol
from stmol import showmol, add_hover
from chembl_structure_pipeline import standardizer


st.title("Chemical Shift Prediction")
st.markdown("#### The chemical shift prediction is done in nmrdb")

# Generate the molecule including drawing option
smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
smiles = st_ketcher(smiles)

mol = Chem.MolFromSmiles(smiles)
mol = Chem.MolToMolBlock(mol)
std_molblock = standardizer.standardize_molblock(mol)


# Add protons is necessary
mol = std_molblock
mol_for_H = Chem.MolFromMolBlock(mol)
mol_with_H = Chem.AddHs(mol_for_H)
smiles_with_H = Chem.MolToSmiles(mol_with_H)

col13C, col1H = st.columns(2)

with col13C:
    create_html_ref_13C = f"https://www.nmrdb.org/service.php?name=nmr-13c-prediction&smiles={smiles}"
    st.page_link(create_html_ref_13C, label="13C Chemical shift prediction (nmrdb)", icon="🤖", help=None, disabled=False, use_container_width=None)
with col1H:
    create_html_ref_1H = f"https://www.nmrdb.org/service.php?name=nmr-1h-prediction&smiles={smiles_with_H}"
    st.page_link(create_html_ref_1H, label="1H Chemical shift prediction (nmrdb)", icon="🤖", help=None, disabled=False, use_container_width=None)
