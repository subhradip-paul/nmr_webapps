import streamlit as st
import pandas as pd
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from streamlit import session_state

rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
from scipy.spatial.transform import Rotation as Rot
from io import StringIO
import py3Dmol
from stmol import showmol, add_hover

script_dir = os.path.dirname(__file__)
csv_file = os.path.join(script_dir, '../dep/NMR_freq_table.csv')

nuctable=pd.read_csv(csv_file)
gyr_ratio_MHz_T=nuctable["GyrHz"]
name_nuc=nuctable["Name"]


def render3dmol(mol, style):
    # Convert RDKit molecule to MOL block format
    mol_block = Chem.MolToMolBlock(mol)

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
    showmol(xyzview,height=500,width=800)

def angle_between_vectors(vec1 , vec2) :
    """
    The function calculates the Euler angles between two vectors
    in the zyz convention. The conventions can be found in scipy's transform module
    :param vec1: a vector of the format [i1, j1, k1]
    :param vec2: b vector of the format [i2, j2, k2]
    :return: a vector of three Euler angles connecting the two vectors
    """
    vec1 = vec1 / np.linalg.norm ( vec1 )
    vec2 = vec2 / np.linalg.norm ( vec2 )
    cp = np.cross ( vec1 , vec2 )
    dp = np.dot ( vec1 , vec2 )
    cosine_ang = dp
    ssc = np.array ( [ [ 0 , -cp[ 2 ] , cp[ 1 ] ] , [ cp[ 2 ] , 0 , -cp[ 0 ] ] , [ -cp[ 1 ] , cp[ 0 ] , 0 ] ] )
    rot_mat = np.identity ( 3 ) + ssc + np.dot ( ssc , ssc ) / (1 + cosine_ang)
    ff = Rot.from_matrix ( rot_mat )
    euler_angles = ff.as_euler ( 'zyz' , degrees='True' )
    return euler_angles




def choose_nucleus(symbol, key_gen):
    """Returns the user's choice when multiple nucleus options exist."""
    matches = nuctable[nuctable["Symbol"] == symbol]

    if matches.shape[0] == 1:
        return matches.iloc[0]["GyrHz"]  # Use the only available match

    elif matches.shape[0] > 1:
        # Let the user choose via a dropdown menu in Streamlit
        choice = st.selectbox(f"Multiple isotopes found for {symbol}. Choose one:",
                              matches["Name"].tolist(), key=str(key_gen) + symbol)
        return matches[matches["Name"] == choice]["GyrHz"].values[0]

    else:
        st.warning(f"No match found for {symbol}. Assigning default value 0.0")
        return 0.0  # Prevents IndexError
def get_gyration_inputs(atom_list):
    """Create UI for isotope selection and return gyration array."""
    gyr_atom = []

    st.subheader("Isotope Selection")

    for idx, symbol in enumerate(atom_list):
        matches = nuctable[nuctable["Symbol"] == symbol]

        if matches.shape[0] == 1:
            gyr_atom.append(matches.iloc[0]["GyrHz"])

        elif matches.shape[0] > 1:
            key = f"isotope_{idx}_{symbol}"

            if key not in st.session_state:
                st.session_state[key] = matches["Name"].iloc[0]

            choice = st.selectbox(
                f"{symbol} (atom {idx+1})",
                matches["Name"].tolist(),
                key=key
            )

            gyr = matches[matches["Name"] == choice]["GyrHz"].values[0]
            gyr_atom.append(gyr)

        else:
            st.warning(f"{symbol} not found. Using 0.")
            gyr_atom.append(0.0)

    return np.array(gyr_atom)

def xyz_to_dipolar_data(xyz_dataframe, gyr_atom):
    num_atoms = len(xyz_dataframe)

    nuc = xyz_dataframe['atom'].tolist()
    coord_xyz = xyz_dataframe[['x', 'y', 'z']].astype(float).to_numpy()

    pl = 6.62607015e-34

    results = []

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):

            r = np.linalg.norm(coord_xyz[i] - coord_xyz[j])

            gyr1 = gyr_atom[i] * 1e6
            gyr2 = gyr_atom[j] * 1e6

            dip = -1e-7 * (gyr1 * gyr2 * pl) / ((r * 1e-10) ** 3)

            euler = np.round(angle_between_vectors(coord_xyz[i], coord_xyz[j]), 2)

            results.append([
                i + 1, nuc[i],
                j + 1, nuc[j],
                round(r, 2),
                round(dip, 2),
                *euler
            ])

    return pd.DataFrame(results, columns=[
        'i', 'Nuc i', 'j', 'Nuc j',
        'Distance (Å)', 'Dipolar Coupling (Hz)',
        'alpha', 'beta', 'gamma'
    ])


def dist2dipole(choice_nuc1, choice_nuc2, distance):
    pl = 6.62607015e-34
    nuc1idx = name_nuc[name_nuc.str.match(choice_nuc1)].index
    nuc2idx = name_nuc[name_nuc.str.match(choice_nuc2)].index

    gyr1=gyr_ratio_MHz_T[nuc1idx[0]]*1e6
    gyr2=gyr_ratio_MHz_T[nuc2idx[0]]*1e6

    dip=-1e-7*(gyr1*gyr2*pl)/((distance*1e-10) ** 3)
    return dip
    

def dipole2dist(choice_nuc1, choice_nuc2, dipole):
    pl = 6.62607015e-34
    nuc1idx = name_nuc[name_nuc.str.match(choice_nuc1)].index
    nuc2idx = name_nuc[name_nuc.str.match(choice_nuc2)].index

    gyr1=gyr_ratio_MHz_T[nuc1idx[0]]*1e6
    gyr2=gyr_ratio_MHz_T[nuc2idx[0]]*1e6

    dist=1e10 * ((1e-7*abs((gyr1*gyr2*pl))/dipole) ** (1/3))
    return dist

def simpson_display(df):
    if 'event' not in st.session_state:
        event = st.dataframe(
            df,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="multi-row",
        )

        st.subheader("For SIMPSON")
        dipole_pairs = event.selection.rows
        filtered_df = df.iloc[dipole_pairs]
        st.code(filtered_df)

def click_button():
    st.session_state.clicked = True

st.title(r"Dipolar Coupling Calculator")
st.divider()

choice_of_calculation = st.radio(":rainbow[**What do you want to calculate?**]", [r"Dipole :arrows_counterclockwise: Distance", r'Dipolar Couplings from Structure'], index=None)

if choice_of_calculation == "Dipole :arrows_counterclockwise: Distance":
    choice_nuc1 = st.selectbox('Nucleus 1:', nuctable.Name.values.tolist(), index=2)
    choice_nuc2 = st.selectbox('Nucleus 2:', nuctable.Name.values.tolist(), index=11)
    choice_calc = st.radio('Method', ('Distance to Dipole', 'Dipole to Distance'))

    if choice_calc == 'Distance to Dipole':
        distance=st.number_input('Enter distance in Angstrom: ', min_value=0.5, value=1.07)
        coup = dist2dipole(choice_nuc1,choice_nuc2,distance)
        st.divider()
        st.write("The dipolar coupling = " + str(np.abs(np.round(coup,2))) + ' Hz')
        st.write("The dipolar coupling = " + str(np.abs(np.round(coup/1e3,2))) + ' kHz')
    else:
        dipole=st.number_input('Enter dipolar coupling in Hz: ', min_value=0.5)
        d=dipole2dist(choice_nuc1,choice_nuc2,dipole)
        st.divider()
        st.latex(r"\text{Distance  = }" + str(np.abs(np.round(d,2))) + r'\ {\AA}')
elif choice_of_calculation == 'Dipolar Couplings from Structure':
    options = st.selectbox("Draw a structure or manually enter *xyz* coordinates", options=["Draw", "Manual"], index=0)
    if options == "Manual":
        xyz_lines = st.file_uploader("Upload an xyz coordinates file", type="xyz")
        if xyz_lines is not None:
            df = pd.read_csv ( xyz_lines , skiprows=2,
                               sep=r"\s+" ,
                               names=[ "atom" , "x" , "y" , "z" ] ,
                               dtype={"atom" : str , "x" : float , "y" : float , "z" : float} )
            st.data_editor(df)
            num_atoms = len(df)
            df_xyz_dipole = xyz_to_dipolar_data ( xyz_dataframe=df , num_atoms=num_atoms )
            st.write ( df_xyz_dipole )

    else:
        smiles = 'cc'
        smiles = st_ketcher(smiles)
        if smiles is not None:
            mol = Chem.MolFromSmiles( smiles)
            if mol is not None:

                # xyz_file = os.path.join(script_dir, '../dep/temp.xyz')
                # Chem.MolToXYZFile(mol, xyz_file)
                remove_1H = st.checkbox('Remove 1H')
                if remove_1H:
                    mol3d = Chem.AddHs(mol)
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 0xf00d  # optional random seed for reproducibility
                    AllChem.EmbedMolecule(mol3d, params)
                    AllChem.MMFFOptimizeMolecule(mol3d)
                    mol3d = Chem.RemoveAllHs(mol3d)
                    mol_to_xyz = mol3d
                    # Dropdown for style selection
                    style = st.selectbox("Select style:", ['stick', 'ball and stick', 'sphere'])
                    render3dmol(mol3d, style)

                else:
                    mol3d = Chem.AddHs(mol)
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 0xf00d  # optional random seed for reproducibility
                    AllChem.EmbedMolecule(mol3d, params)
                    AllChem.MMFFOptimizeMolecule(mol3d)
                    mol_to_xyz = mol3d
                    # Dropdown for style selection
                    style = st.selectbox("Select style:", ['stick', 'ball and stick', 'sphere'])
                    render3dmol(mol3d, style)


                xyz_string = Chem.MolToXYZBlock(mol_to_xyz)
                xyz_lines = xyz_string.strip().splitlines()  # Split into lines and remove empty spaces
                num_atoms = int(xyz_lines[0])

            if "run_calc" not in st.session_state:
                st.session_state.run_calc = False

            if st.button("Calculate"):
                st.session_state.run_calc = True

            if st.session_state.run_calc:
                df = pd.read_csv(
                    StringIO("\n".join(xyz_lines[2:])),
                    sep=r"\s+",
                    names=["atom", "x", "y", "z"]
                )

                st.write("### Input structure")
                st.dataframe(df)

                # 👉 UI phase (stable)
                gyr_atom = get_gyration_inputs(df["atom"].tolist())

                # 👉 Compute phase
                df_xyz_dipole = xyz_to_dipolar_data(df, gyr_atom)

                st.write("### Dipolar couplings")
                st.dataframe(df_xyz_dipole, hide_index=True)







