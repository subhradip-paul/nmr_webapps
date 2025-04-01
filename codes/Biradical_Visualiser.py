import os

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
import rdkit.Chem.Descriptors as MolDescriptors
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
from stmol import showmol, makeobj, add_hover
from chembl_structure_pipeline import standardizer
import py3Dmol



def return_smiles_code_mol():
    molecule_dict = {
        "bTbK": "CC1(C)CC2(CC(C)(C)N1[O])OCC1(CO2)COC2(CC(C)(C)N([O])C(C)(C)C2)OC1",
        "bTurea": "CC1(C)CC(NC(=O)NC2CC(C)(C)N([O])C(C)(C)C2)CC(C)(C)N1[O]",
        "bCTbK": "[O]N1C2(CCCCC2)CC2(CC13CCCCC3)OCC1(CO2)COC2(CC3(CCCCC3)N([O])C3(CCCCC3)C2)OC1",
        "bMTbK": "CC1CC(C)CC2(C1)CC1(CC3(CC(C)CC(C)C3)N2[O])OCC2(CO1)COC1(CC3(CC(C)CC(C)C3)N([O])C3(CC(C)CC(C)C3)C1)OC2",
        "ABK": "[O]N1C2CCCC1CC1(C2)OCC2(CO1)COC1(CC3CCCC(C1)N3[O])OC2",
        "ABU": "[O]N1C2CCCC1CC(NC(=O)NC1CC3CCCC(C1)N3[O])C2",
        "TEKPol": "[O]N1C2(CCC(C3=CC=CC=C3)CC2)CC2(CC13CCC(C1=CC=CC=C1)CC3)OCC1(CO2)COC2(CC3(CCC(C4=CC=CC=C4)CC3)N([O])C3(CCC(C4=CC=CC=C4)CC3)C2)OC1",
        "PyPol": "[2H]N(C(=O)NC1CC2(CCOCC2)N([O])C2(CCOCC2)C1)C1CC2(CCOCC2)N([O])C2(CCOCC2)C1",
        "AMUPol": "[H]N(C(=O)N(CCOCCOCCOCCOC)C1CC2(CCOCC2)N([O])C2(CCOCC2)C1)C1CC2(CCOCC2)N([O])C2(CCOCC2)C1",
        "TEKPolCbm": "COC(=O)N(C(=O)OC1CC2(CCC(C3=CC=CC=C3)CC2)N([O])C2(CCC(C3=CC=CC=C3)CC2)C1)C1CC2(CCC(C3=CC=CC=C3)CC2)N([O])C2(CCC(C3=CC=CC=C3)CC2)C1",
        "PyPolCbm": "[H]N(C(=O)OC1CC2(CCOCC2)N([O])C2(CCOCC2)C1)C1CC2(CCOCC2)N([O])C2(CCOCC2)C1",
        "AMUPolCbm": "N(CCOCOCOCOC)(C1CC2(CCOCC2)N([O])C2(CCOCC2)C1)C(OC1CC2(CCOCC2)N([O])C2(CCOCC2)C1)=O",
        "AsymPolPOK": "[K+][O-]P(=O)([O-][K+])OC1CCC2(CC(NC(=O)C3=CC(C)(C)N([O])C3(C)C)CC3(CCC(OP(=O)([O-][K+])[O-][K+])CC3)N2[O])CC1 |^1:22,41|",
        "AsymPolTEK": "CC1(C)C=C(C(=O)NC2CC3(CCC(C4=CC=CC=C4)CC3)N([O])C3(CCC(C4=CC=CC=C4)CC3)C2)C(C)(C)N1[O]",
        "HydroPol": "CC[11CH2]N(C(=O)N(C[C](C)O)C1C[C@]2(CC(C)OC(C)C2)N([O])[C@@]2(CC(C)OC(C)C2)C1)C1C[C@]2(CC(C)OC(C)C2)N([O])[C@@]2(CC(C)OC(C)C2)C1",
        "TEKPolCbo": "[O]N1C2(CCC(C3=CC=CC4=C3C=CC=C4)CC2)CC2(CC13CCC(C1=CC=CC4=C1C=CC=C4)CC3)OCC1(CO2)COC2(CC3(CCC(C4=CC=CC5=C4C=CC=C5)CC3)N([O])C3(CCC(C4=CC=CC5=C4C=CC=C5)CC3)C2)OC1",
        "NaphPol": "[O]N1C2(CCC(C3=CC=CC=C3)CC2)CC(OC(=O)OC2CC3(CCC(C4=CC=CC=C4)CC3)N([O])C3(CCC(C4=CC=CC=C4)CC3)C2)CC12CCC(C1=CC=CC=C1)CC2",
        "C-bcTol": "O=C(NC1CC2(CCC(O)CC2)N([O])C2(CCC(O)CC2)C1)NC1CC2(CCC(O)CC2)N([O])C2(CCC(O)CC2)C1",
        "4-Oxo-TEMPO": "CC1(C)CC(=O)CC(C)(C)N1[O]",
        "TEMPO": "CC1(C)CCCC(C)(C)N1[O]",
        "Trityl": "C1=CC=C(C=C1)[C](C1=CC=CC=C1)C1=CC=CC=C1",
        "TOTAPOL": "CC1(C)CC(CC(C)(C)N1[O])NCC(COC2CC(C)(C)N(C(C)(C)C2)[O])O",
        "cAsymPolTEK": "C1C=CC(C2CCC3(N([O])C4(CCC(C5C=CC=CC=5)CC4)CC(NC(C4C5(CCCCC5)N([O])C5(CCCCC5)C=4)=O)C3)CC2)=CC=1",
        "Ox063": "[O-]C(=O)C1=C2C(SC%87%88S2)=C([C@](C2=C3C(SC%89%90S3)=C(C(=O)[O-])C3=C2SC%91%92S3)C2=C3C(SC%93%94S3)=C(C(=O)[O-])C3=C2SC%95%96S3)C2=C1SC%97%98S2.[*:1]%97.[*:1]%98.[*:1]%95.[*:1]%96.[*:1]%93.[*:1]%94.[*:1]%91.[*:1]%92.[*:1]%89.[*:1]%90.[*:1]%87.[*:1]%88 |^1:10,$;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;_R1;_R1;_R1;_R1;_R1;_R1;_R1;_R1;_R1;_R1;_R1;_R1$|",
        "PyrrhoTriPol": "[O-]C(=O)C1=C2C(SC(C)(C)S2)=C([C@](C2=C3C(SC(C)(C)S3)=C(C(=O)[O-])C3=C2SC(C)(C)S3)C2=C3C(SC(C)(C)S3)=C(C(=O)ON3CCN(CC3)C(=O)C3C(N([O])C(C=3)(C)C)(C)C)C3=C2SC(C)(C)S3)C2=C1SC(C)(C)S2 |^1:12,55|",
        "PyrrhoTriPol-OMe": "CC1(Sc2c(c3SC(C)(C)Sc3c([C](c3c4SC(Sc4c(C(ON4CCN(C(C5C(C)(C)N([O])C(C)(C)C=5)=O)CC4)=O)c4SC(Sc34)(C)C)(C)C)c3c4SC(Sc4c(C(OC)=O)c4SC(Sc34)(C)C)(C)C)c2S1)C(OC)=O)C |^1:13,33|",
        "BDPA": "[C](/C(/C1C=CC=CC=1)=C1\C2C=CC=CC=2C2C=CC=CC\\1=2)1C2C=CC=CC=2C2C=CC=CC1=2 |^1:0|",
    }

    # Convert dictionary to DataFrame
    df = pd.DataFrame(list(molecule_dict.items()), columns=["Biradical Name", "SMILES Code"])
    df.sort_values(by = "Biradical Name", inplace = True, key = lambda col: col.str.lower())
    return df

def plot_2d_struct_span_ketcher(smiles_2dfn):
    st.subheader("You can edit the structure of the molecule here")
    st.write('Credit: [Streamlit-Ketcher](https://github.com/mik-laj/streamlit-ketcher)')
    ketcher_window_smiles = st_ketcher(smiles_2dfn)
    mol = Chem.MolFromSmiles( ketcher_window_smiles )
    mol = Chem.MolToMolBlock(mol)
    std_molblock = standardizer.standardize_molblock(mol)
    mol = Chem.MolFromMolBlock(std_molblock)
    mol = Chem.RemoveHs(mol)
    smiles_modified = Chem.MolToSmiles(mol)

    st.subheader("You can see the structure of the molecule here")
    if mol:
        drawing_options = Draw.MolDrawOptions()
        drawing_options.addStereoAnnotation = True
        drawing_options.includeAtomTags = True
        drawing_options.addAtomIndices = True
        img = Draw.MolToImage(mol,size=(400,400), options=drawing_options)
        st.image(img, caption=f"{biradical_chosen}", use_container_width='auto')
        return mol, smiles_modified
    else:
        raise Exception("Not a valid molecule")



def generate_3d_molecule(smiles_3dfn, style):
    try:
        mol = Chem.MolFromSmiles( smiles_3dfn )
        mol3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3 ()
        params.randomSeed = 0xf00d  # optional random seed for reproducibility
        AllChem.EmbedMolecule ( mol3d , params )
        AllChem.MMFFOptimizeMolecule ( mol3d )

        # Convert to MolBlock
        mol_no_h = Chem.RemoveAllHs(mol3d)
        mol_block = Chem.MolToMolBlock(mol_no_h)
        xyzview = py3Dmol.view()
        xyzview.addModel(mol_block, 'mol')  # Pass MOL block instead of RDKit Mol object

        # Set molecular visualization style
        if style == 'stick':
            xyzview.setStyle({'stick': {}})
        elif style == 'ball and stick':
            xyzview.setStyle({'sphere': {'scale': 0.3}, 'stick': {}})
        elif style == 'sphere':
            xyzview.setStyle({'sphere': {}})

        colour_bg = st.color_picker('Pick a background colour', '#9DEBF7')
        xyzview.setBackgroundColor(colour_bg)
        xyzview.zoomTo()
        add_hover(xyzview)
        showmol(xyzview, height=500, width=800)

        mol_xyz_noh = Chem.MolToXYZBlock(mol_no_h)
        mol_xyz = Chem.MolToXYZBlock(mol3d)

        return mol_xyz_noh, mol_xyz
        
    except Exception as e:
        st.error(f"Error generating 3D molecule: {e}")
        return None

# Load SMILES data

df_smile_file = return_smiles_code_mol()
biradical_choice = df_smile_file["Biradical Name"]

# Streamlit app
st.title("Biradical Structure Viewer")

# Select biradical
biradical_chosen = st.selectbox("Choose a biradical", options=biradical_choice)
smiles = df_smile_file.loc[df_smile_file["Biradical Name"] == biradical_chosen, "SMILES Code"].values[0]

st.divider()
st.header("2D Structure Viewer")
mol_2d, ketcher_modified_smiles = plot_2d_struct_span_ketcher(smiles)


# 3D structure viewer
st.divider()
st.header("3D Structure Viewer")
st.markdown('''
This is a 3D structure generated from the .mol file by RDKit.
It is not a DFT optimised structure. Please consider it as a representative structure
The 3D view is generated using the package stmol:
Nápoles-Duarte JM, Biswas A, Parker MI, Palomares-Baez JP, Chávez-Rojo MA and Rodríguez-Valdez LM (2022) Stmol: A component for building interactive molecular visualizations within streamlit web-applications. Front. Mol. Biosci. 9:990846. doi: 10.3389/fmolb.2022.990846
''')
style = st.selectbox("Select style:", ['stick', 'ball and stick', 'sphere'])
xyz_without_H, xyz = generate_3d_molecule(ketcher_modified_smiles, style)

# Creating options to download the xyz coordinates

col_xyz, col_xyz_noH = st.columns(2)

with col_xyz_noH:
    if st.button("Show XYZ coordinates without H"):
        st.code(xyz_without_H)
        st.download_button('Download the coordinates as a .xyz file', xyz_without_H, file_name=f'{biradical_chosen}_withoutH.xyz', mime='.xyz', key='biradical without H.xyz')
with col_xyz:
    if st.button("Show XYZ coordinates with H"):
        st.code(xyz)
        st.download_button('Download the coordinates as a .xyz file', xyz, file_name=f'{biradical_chosen}_withH.xyz', mime='.xyz', key='biradical.xyz')

st.divider()
st.header("Some Important Properties of the Biradical Structure")
molecular_weight = MolDescriptors.ExactMolWt(mol_2d)
num_rotable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_2d, strict=False)

num_radicals = MolDescriptors.NumRadicalElectrons(mol_2d)

st.subheader(f'{biradical_chosen}')
col1, col2 = st.columns(2)

with col1:
    st.metric(label="Molecular Weight", value=round(molecular_weight, 2))
with col2:
    st.metric(label="Number of Rotatable Bonds", value=num_rotable_bonds)

st.metric(label="Number of Radical Electrons", value=num_radicals)

# Calculate all molecular descriptors
mol_props = MolDescriptors.CalcMolDescriptors(mol_2d)

# Convert dictionary keys and values into a list of tuples for display
mol_data = [(key, round(value, 3) if isinstance(value, (int, float)) else value) for key, value in mol_props.items()]

# Streamlit UI
st.markdown("##### Rest of the Molecular Properties")

# Display properties in a table
st.dataframe(mol_data, column_config={"0": "Property", "1": "Value"}, hide_index=True, use_container_width=True)

