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

def xyz_file_to_dipolar_data(xyzfile) :
    """
    The function takes an .xyz file of a molecular structure and calculates the
    dipolar coupling and Euler angle between the different nuclei in the principal axis frame.
    :param xyzfile: Molecular structure file of the format .xyz
    :return: a pandas Dataframe containing the pair of nuclei, the dipolar coupling in Hz,
    and the Euler angles between the two tensors.
    """

    mol = pd.read_csv ( xyzfile , sep=r'\s+' , skiprows=2 , names=[ 'atom' , 'x' , 'y' , 'z' ] , index_col=False )
    nuc = mol['atom'].tolist()
    coord_xyz = mol[ [ 'x' , 'y' , 'z' ] ].to_numpy ()
    pl = 6.62607015e-34
    df_xyz_to_dip = pd.DataFrame ( columns=[ 'i' , 'j' , 'dip' , 'alpha' , 'beta' , 'gamma' ] )

    gyr_atom = np.zeros(len(nuc))

    for idx, nucleus in enumerate(nuc):
        gyr_atom[idx] = nuctable[nuctable.Name.isin([nucleus]) == True]['GyrHz'].values



    dist = np.zeros ( [ np.shape ( coord_xyz )[ 0 ] , np.shape ( coord_xyz )[ 0 ] ] )
    dip = np.zeros ( [ np.shape ( coord_xyz )[ 0 ] , np.shape ( coord_xyz )[ 0 ] ] )
    for idx in range ( 0 , np.shape ( coord_xyz )[ 0 ] ) :
        gyr1 = gyr_atom[ idx ] * 1e6
        for j in range ( idx + 1 , np.shape ( coord_xyz )[ 0 ] ) :
            dist[ idx ][ j ] = np.sqrt ( np.sum ( (coord_xyz[ idx ] - coord_xyz[ j ]) ** 2 ) )
            gyr2 = gyr_atom[ j ] * 1e6
            dip[ idx ][ j ] = -1e-7 * (gyr1 * gyr2 * pl) / ((dist[ idx ][ j ] * 1e-10) ** 3)

    for idx in range ( 0 , np.shape ( coord_xyz )[ 0 ] ) :
        for j in range ( idx + 1 , np.shape ( coord_xyz )[ 0 ] ) :
            euler_angles = np.round ( angle_between_vectors ( coord_xyz[ idx ] , coord_xyz[ j ] ) , 2 )
            df_xyz_to_dip.loc[ idx * np.shape ( coord_xyz )[ 0 ] + j ] = [ idx + 1 , j + 1 ,
                                                                           np.round ( dip[ idx ][ j ] , 2 ) ,
                                                                           *euler_angles ]

    return df_xyz_to_dip



def dist2dipole(choice_nuc1, choice_nuc2, distance):
    pl = 6.62607015e-34
    nuc1idx = name_nuc[name_nuc.str.match(choice_nuc1)].index
    nuc2idx = name_nuc[name_nuc.str.match(choice_nuc2)].index

    gyr1=gyr_ratio_MHz_T[nuc1idx[0]]*1e6
    gyr2=gyr_ratio_MHz_T[nuc2idx[0]]*1e6

    dip=-1e-7*(gyr1*gyr2*pl)/((distance*1e-10) ** 3);
    return dip
    # print("The dipolar coupling = " + str(np.abs(np.round(dip,2))) + ' Hz')
    # print("The dipolar coupling = " + str(np.abs(np.round(dip/1e3,2))) + ' kHz')
    

def dipole2dist(choice_nuc1, choice_nuc2, dipole):
    pl = 6.62607015e-34
    nuc1idx = name_nuc[name_nuc.str.match(choice_nuc1)].index
    nuc2idx = name_nuc[name_nuc.str.match(choice_nuc2)].index

    gyr1=gyr_ratio_MHz_T[nuc1idx[0]]*1e6
    gyr2=gyr_ratio_MHz_T[nuc2idx[0]]*1e6

    dist=1e10 * ((1e-7*abs((gyr1*gyr2*pl))/dipole) ** (1/3))
    return dist
    # print("The distance  = " + str(np.abs(np.round(dist,2))) + ' A')
    

st.title(r"Dipole :arrows_counterclockwise: Distance")
st.divider()

choice_of_calculation = st.radio("What do you want to calculate?", [r"Dipole :arrows_counterclockwise: Distance", r'Dipolar Couplings from Structure :molecule:'], index=None)

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
elif choice_of_calculation == 'Dipolar Couplings from Structure :molecule:':
    smiles = 'cc'
    ketcher_window_smiles = st_ketcher(smiles)
    mol = Chem.MolFromSmiles( ketcher_window_smiles)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d  # optional random seed for reproducibility
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    xyz_file = os.path.join(script_dir, '../dep/temp.xyz')
    Chem.MolToXYZFile(mol, xyz_file)

