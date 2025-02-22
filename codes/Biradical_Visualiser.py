import os

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
import rdkit.Chem.Descriptors as MolDescriptors
rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
from stmol import showmol, makeobj, add_hover


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
        "AsymPolPOK": "[O-]P(OC1CCC2(N([O])C3(CCC(OP([O-])([O-])=O)CC3)CC(NC(C3C(C)(C)N([O])C(C)(C)C=3)=O)C2)CC1)([O-])=O",
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
    }

    # Convert dictionary to DataFrame
    df = pd.DataFrame(list(molecule_dict.items()), columns=["Biradical Name", "SMILES Code"])
    df.sort_values(by = "Biradical Name", inplace = True, key = lambda col: col.str.lower())
    return df

def plot_2d_struct_span_ketcher(smiles_2dfn):
    st.subheader("You can edit the structure of the molecule here")
    st.write('Credit: [Streamlit-Ketcher](https://github.com/mik-laj/streamlit-ketcher)')
    mol = Chem.MolFromSmiles(smiles_2dfn)
    mol = Chem.RemoveHs(mol)
    smiles_noh = Chem.MolToSmiles(mol)
    ketcher_window_smiles = st_ketcher(smiles_noh)
    mol = Chem.MolFromSmiles( ketcher_window_smiles )
    st.subheader("You can see the structure of the molecule here")
    if mol:
        img = Draw.MolToImage(mol,size=(400,400))
        st.image(img, caption=f"{biradical_chosen}", use_container_width=True)
        return mol, ketcher_window_smiles
    else:
        raise Exception("Not a valid molecule")



def generate_3d_molecule(smiles_3dfn):
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
        py3dmol_obj_return = makeobj ( mol_block , molformat='mol' , style='stick' , background='white' )
        add_hover ( py3dmol_obj_return , backgroundColor='white' , fontColor='black' )
        return py3dmol_obj_return
        
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

py3dmol_obj = generate_3d_molecule(ketcher_modified_smiles)
showmol(py3dmol_obj, height=500, width=800)
st.divider()
st.header("Some Important Properties of the Biradical Structure")
molecular_weight = MolDescriptors.ExactMolWt(mol_2d)
num_aromatic_rings = MolDescriptors.NumAromaticRings(mol_2d)
num_aliphatic_rings = MolDescriptors.NumRadicalElectrons(mol_2d)
num_rings = MolDescriptors.NumAliphaticRings(mol_2d)
num_radicals = MolDescriptors.NumRadicalElectrons(mol_2d)

st.subheader(f'{biradical_chosen}')
col1, col2 = st.columns(2)

with col1:
    st.metric(label="Molecular Weight", value=round(molecular_weight, 2))
    st.metric(label="Number of Aromatic Rings", value=num_aromatic_rings)

with col2:
    st.metric(label="Number of Aliphatic Rings", value=num_aliphatic_rings)
    st.metric(label="Total Number of Rings", value=num_rings)

st.metric(label="Number of Radical Electrons", value=num_radicals)

# Calculate all molecular descriptors
mol_props = MolDescriptors.CalcMolDescriptors(mol_2d)

# Convert dictionary keys and values into a list of tuples for display
mol_data = [(key, round(value, 3) if isinstance(value, (int, float)) else value) for key, value in mol_props.items()]

# Streamlit UI
st.markdown("##### Rest of the Molecular Properties")

# Display properties in a table
st.dataframe(mol_data, column_config={"0": "Property", "1": "Value"}, hide_index=True, use_container_width=True)

