import streamlit as st
from rdkit import Chem
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
import polars as pl
from chembl_structure_pipeline import standardizer
from rdkit.Chem import Draw

st.title("Chemical Shift")

st.header("Isotropic Chemical Shift Prediction")

st.markdown("#### The chemical shift prediction is done in NMRium")

# Generate the molecule including drawing option
smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'

smiles_from_ketcher = st_ketcher(smiles)
if smiles_from_ketcher is not None:
    smiles = smiles_from_ketcher

mol = Chem.MolFromSmiles(smiles)
mol = Chem.MolToMolBlock(mol)
std_molblock = standardizer.standardize_molblock(mol)


# Add protons is necessary
mol = std_molblock
mol_for_H = Chem.MolFromMolBlock(mol)
mol_with_H = Chem.AddHs(mol_for_H)
smiles_with_H = Chem.MolToSmiles(mol_with_H)

st.info('Before chemical shift prediction press the apply button')



# NMRium 13C and 1H prediction
nmrium_13C_url = f"https://www.nmrium.org/predict?smiles={smiles_with_H}"
st.link_button("Chemical shift prediction (NMRium)", nmrium_13C_url,
               help='Clicking here will open an external page where the chemical shifts will be predicted',
               type='primary', icon=":material/arrow_circle_down:", disabled=False, use_container_width=True)



st.warning('The chemical shifts are mostly for solution NMR')

st.info('For other nuclei and for more comprehensive list of molecules, visit the website: https://nmrshiftdb.nmr.uni-koeln.de')


st.header("Chemical Shift Anisotropy")

st.write('The chemical shift data is taken from the publication: Chen, Xi, and Chang-Guo Zhan. '
         '“First-Principles Studies of C-13 NMR Chemical Shift Tensors of Amino Acids in Crystal State.” '
         'Journal of Molecular Structure: THEOCHEM 682, (2004): 73–82.[Link](https://doi.org/10.1016/j.theochem.2004.05.027.)')


# Data for the 20 standard amino acids
data = [
    ["Alanine", "Ala", "A", "CC(C(=O)O)N"],
    ["Arginine", "Arg", "R", "NC(CCCN=C(N)N)C(=O)O"],
    ["Asparagine", "Asn", "N", "NC(CC(=O)N)C(=O)O"],
    ["Aspartic acid", "Asp", "D", "NC(CC(=O)O)C(=O)O"],
    ["Cysteine", "Cys", "C", "NC(CS)C(=O)O"],
    ["Glutamic acid", "Glu", "E", "NC(CCC(=O)O)C(=O)O"],
    ["Glutamine", "Gln", "Q", "NC(CCC(=O)N)C(=O)O"],
    ["Glycine", "Gly", "G", "NCC(=O)O"],
    ["Histidine", "His", "H", "NC(CC1=CN=CN1)C(=O)O"],
    ["Isoleucine", "Ile", "I", "CC[C@H](C)[C@H](N)C(=O)O"],
    ["Leucine", "Leu", "L", "CC(C)CC(N)C(=O)O"],
    ["Lysine", "Lys", "K", "NC(CCCCN)C(=O)O"],
    ["Methionine", "Met", "M", "CSCC(N)C(=O)O"],
    ["Phenylalanine", "Phe", "F", "NC(CC1=CC=CC=C1)C(=O)O"],
    ["Proline", "Pro", "P", "OC1CCCN1"],
    ["Serine", "Ser", "S", "NC(CO)C(=O)O"],
    ["Threonine", "Thr", "T", "CC(O)C(N)C(=O)O"],
    ["Tryptophan", "Trp", "W", "NC(CC1=CNC2=CC=CC=C12)C(=O)O"],
    ["Tyrosine", "Tyr", "Y", "NC(CC1=CC=C(O)C=C1)C(=O)O"],
    ["Valine", "Val", "V", "CC(C)C(N)C(=O)O"],
]

# Create Polars DataFrame
df = pl.DataFrame(
    data,
    schema=["Amino Acid", "3-Letter Code", "1-Letter Code", "SMILES"], orient = 'row',
)

amino_acid_selection = st.selectbox("Choose an amino acid", options=df["Amino Acid"], index=7)

if amino_acid_selection is not None:
    mask = df.filter(pl.col("Amino Acid")==amino_acid_selection)
    mol = Chem.MolFromSmiles(mask['SMILES'].item())
    if mol:
        drawing_options = Draw.MolDrawOptions()
        drawing_options.addStereoAnnotation = True
        drawing_options.includeAtomTags = True
        drawing_options.addAtomIndices = False
        img = Draw.MolToImage(mol,size=(400,400), options=drawing_options)
        st.image(img, caption=f"{amino_acid_selection}", use_container_width='auto')

import os
script_dir = os.path.dirname(__file__)
csv_file = os.path.join(script_dir, '../dep/csa_amino_acids.csv')


df_csa_tensor_amino_acid = pl.read_csv(csv_file)

df_csa_tensor_amino_acid_sel=df_csa_tensor_amino_acid.filter(pl.col("Amino Acid").str.contains(amino_acid_selection))



# Compute `diso`
df_csa_tensor_amino_acid_sel = df_csa_tensor_amino_acid_sel.with_columns(
    ((df_csa_tensor_amino_acid_sel["δ11"] + df_csa_tensor_amino_acid_sel["δ22"] + df_csa_tensor_amino_acid_sel["δ33"]) / 3).alias("Trace")
)
# Use `pl.when().then().otherwise()` for conditional logic
df_csa_tensor_amino_acid_sel = df_csa_tensor_amino_acid_sel.with_columns(
    pl.when(abs(df_csa_tensor_amino_acid_sel["δ11"]-df_csa_tensor_amino_acid_sel["Trace"]) >= abs(df_csa_tensor_amino_acid_sel["δ33"]-df_csa_tensor_amino_acid_sel["Trace"]))
    .then((-df_csa_tensor_amino_acid_sel["Trace"] + df_csa_tensor_amino_acid_sel["δ11"]) * 1.5)
    .when(abs(df_csa_tensor_amino_acid_sel["δ11"]-df_csa_tensor_amino_acid_sel["Trace"]) <= abs(df_csa_tensor_amino_acid_sel["δ33"]-df_csa_tensor_amino_acid_sel["Trace"]))
    .then((-df_csa_tensor_amino_acid_sel["Trace"] + df_csa_tensor_amino_acid_sel["δ33"]) * 1.5)
    .otherwise(0)
    .alias("CSA (ppm)"),
    pl.when(abs(df_csa_tensor_amino_acid_sel["δ11"]-df_csa_tensor_amino_acid_sel["Trace"]) >= abs(df_csa_tensor_amino_acid_sel["δ33"]-df_csa_tensor_amino_acid_sel["Trace"]))
    .then((df_csa_tensor_amino_acid_sel["δ22"] - df_csa_tensor_amino_acid_sel["δ33"]) / (df_csa_tensor_amino_acid_sel["δ11"] - df_csa_tensor_amino_acid_sel["Trace"]))
    .otherwise((df_csa_tensor_amino_acid_sel["δ22"] - df_csa_tensor_amino_acid_sel["δ11"]) / (df_csa_tensor_amino_acid_sel["δ33"] - df_csa_tensor_amino_acid_sel["Trace"]))
    .alias("η"),
)



st.dataframe(df_csa_tensor_amino_acid_sel, column_config={
        "Trace": st.column_config.NumberColumn(
            "Tr (δ11 + δ22 + δ33)",
            help="Haeberlen Mehring Convention",
            format="%.2f",
        ),
    "η":st.column_config.NumberColumn(
            "η",
            help="Anisotropy Mehring Convention",
            format="%.2f",
        ),
}
        )

df_selected_amino_acid = df_csa_tensor_amino_acid_sel.filter(df_csa_tensor_amino_acid_sel["Amino Acid"].str.contains(amino_acid_selection) == True)






