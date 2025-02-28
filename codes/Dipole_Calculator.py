import streamlit as st
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
from rdkit.Chem import Draw


script_dir = os.path.dirname(__file__)
csv_file = os.path.join(script_dir, '../dep/NMR_freq_table.csv')

nuctable=pd.read_csv(csv_file)
gyr_ratio_MHz_T=nuctable["GyrHz"]
name_nuc=nuctable["Name"]


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

# def xyz_file_to_dipolar_data(xyz_dataframe, num_atoms) :
#     """
#     The function takes an .xyz file of a molecular structure and calculates the
#     dipolar coupling and Euler angle between the different nuclei in the principal axis frame.
#     :param xyz_dataframe: Molecular structure file of the format of a dataframe
#     :return: a pandas Dataframe containing the pair of nuclei, the dipolar coupling in Hz,
#     and the Euler angles between the two tensors.
#     """
#
#     # mol = pd.read_csv ( xyzfile , sep=r'\s+' , skiprows=2 , names=[ 'atom' , 'x' , 'y' , 'z' ] , index_col=False )
#     mol = xyz_dataframe
#     nuc = mol['atom'].tolist()
#     coord_xyz = mol[ [ 'x' , 'y' , 'z' ] ].astype( np.float64 ).to_numpy()
#     pl = 6.62607015e-34
#     df_xyz_to_dip = pd.DataFrame ( columns=[ 'i' , 'Nuc i',  'j', 'Nuc j' , 'Distance', 'Dipolar Coupling' , 'alpha' , 'beta' , 'gamma' ] )
#
#     gyr_atom = np.zeros(int(num_atoms))
#     for idx, nucleus in enumerate(nuc):
#         # gyr_atom[idx] = nuctable[nuctable.Symbol.isin([nucleus]) == True]['GyrHz'].values
#         gyr_atom[idx] = choose_nucleus(nucleus)
#
#
#
#
#     dist = np.zeros ( [ np.shape ( coord_xyz )[ 0 ] , np.shape ( coord_xyz )[ 0 ] ] )
#     dip = np.zeros ( [ np.shape ( coord_xyz )[ 0 ] , np.shape ( coord_xyz )[ 0 ] ] )
#     for idx in range ( 0 , np.shape ( coord_xyz )[ 0 ] ) :
#         gyr1 = gyr_atom[ idx ] * 1e6
#         for j in range ( idx + 1 , np.shape ( coord_xyz )[ 0 ] ) :
#             dist[ idx ][ j ] = np.sqrt ( np.sum ( (coord_xyz[ idx ] - coord_xyz[ j ]) ** 2 ) )
#             gyr2 = gyr_atom[ j ] * 1e6
#             dip[ idx ][ j ] = -1e-7 * (gyr1 * gyr2 * pl) / ((dist[ idx ][ j ] * 1e-10) ** 3)
#
#     for idx in range ( 0 , np.shape ( coord_xyz )[ 0 ] ) :
#         for j in range ( idx + 1 , np.shape ( coord_xyz )[ 0 ] ) :
#             euler_angles = np.round ( angle_between_vectors ( coord_xyz[ idx ] , coord_xyz[ j ] ) , 2 )
#             df_xyz_to_dip.loc[ idx * np.shape ( coord_xyz )[ 0 ] + j ] = [ idx , nuc[idx], j , nuc[j], np.round ( dist[idx][j] , 2 ),
#                                                                            np.round ( dip[ idx ][ j ] , 2 ) ,
#                                                                            *euler_angles ]
#
#     return df_xyz_to_dip


def choose_nucleus(symbol):
    """Returns the user's choice when multiple nucleus options exist."""
    matches = nuctable[nuctable["Symbol"] == symbol]

    if matches.shape[0] == 1:
        return matches.iloc[0]["GyrHz"]  # Use the only available match

    elif matches.shape[0] > 1:
        # Let the user choose via a dropdown menu in Streamlit
        choice = st.selectbox(f"Multiple isotopes found for {symbol}. Choose one:",
                              matches["Name"].tolist(), key=symbol)
        return matches[matches["Name"] == choice]["GyrHz"].values[0]

    else:
        st.warning(f"No match found for {symbol}. Assigning default value 0.0")
        return 0.0  # Prevents IndexError


def xyz_file_to_dipolar_data(xyz_dataframe, num_atoms):
    """
    The function takes an .xyz file of a molecular structure and calculates the
    dipolar coupling and Euler angle between different nuclei in the principal axis frame.

    :param xyz_dataframe: DataFrame containing atomic symbols and coordinates
    :param num_atoms: Number of atoms in the structure
    :return: A pandas DataFrame with nuclear pairs, dipolar coupling, and Euler angles
    """

    mol = xyz_dataframe
    nuc = mol['atom'].tolist()
    coord_xyz = mol[['x', 'y', 'z']].astype(np.float64).to_numpy()
    pl = 6.62607015e-34  # Planck constant

    # DataFrame to store dipolar couplings and Euler angles
    df_xyz_to_dip = pd.DataFrame(columns=['i', 'Nuc i', 'j', 'Nuc j', 'Distance (Å)',
                                          'Dipolar Coupling (Hz)', 'alpha', 'beta', 'gamma'])

    # Assign gyromagnetic ratios, handling missing nuclei
    gyr_atom = np.zeros(int(num_atoms))
    for idx, nucleus in enumerate(nuc):
        gyr_atom[idx] = choose_nucleus(nucleus)
        if gyr_atom[idx] == 0.0:
            st.warning(f"Warning: Nucleus {nucleus} not found in the database. Dipolar interaction may be incorrect.")

    # Initialize distance and dipole matrices
    dist = np.zeros((num_atoms, num_atoms))
    dip = np.zeros((num_atoms, num_atoms))

    # Compute distances and dipolar couplings
    for idx in range(num_atoms):
        gyr1 = gyr_atom[idx] * 1e6  # Convert to Hz
        for j in range(idx + 1, num_atoms):
            dist[idx, j] = np.linalg.norm(coord_xyz[idx] - coord_xyz[j])  # Distance in Å
            gyr2 = gyr_atom[j] * 1e6
            dip[idx, j] = -1e-7 * (gyr1 * gyr2 * pl) / ((dist[idx, j] * 1e-10) ** 3)  # Dipolar coupling in Hz

    # Store results in DataFrame
    for idx in range(num_atoms):
        for j in range(idx + 1, num_atoms):
            euler_angles = np.round(angle_between_vectors(coord_xyz[idx], coord_xyz[j]), 2)
            df_xyz_to_dip.loc[len(df_xyz_to_dip)] = [
                idx + 1, nuc[idx], j + 1, nuc[j], round(dist[idx, j], 2),
                round(dip[idx, j], 2), *euler_angles
            ]

    return df_xyz_to_dip


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
            else:
                mol3d = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 0xf00d  # optional random seed for reproducibility
                AllChem.EmbedMolecule(mol3d, params)
                AllChem.MMFFOptimizeMolecule(mol3d)
                mol_to_xyz = mol3d

            # Draw the final molecule
            drawing_options = Draw.MolDrawOptions()
            drawing_options.addStereoAnnotation = True
            drawing_options.includeAtomTags = True
            img = Draw.MolToImage(mol3d, size=(400, 400), options=drawing_options)
            st.image(img, use_container_width='auto')


            xyz_string = Chem.MolToXYZBlock(mol_to_xyz)
            xyz_lines = xyz_string.strip().splitlines()  # Split into lines and remove empty spaces
            num_atoms = int(xyz_lines[0])

            df = pd.read_csv(StringIO("\n".join(xyz_lines[2:])),
                             delim_whitespace=True,
                             names=["atom", "x", "y", "z"],
                             dtype={"atom": str, "x": float, "y": float, "z": float})
            st.write(df)
            df_xyz_dipole =  xyz_file_to_dipolar_data(xyz_dataframe=df, num_atoms=num_atoms)
            st.write(df_xyz_dipole)




