import streamlit as st
import pandas as pd
import numpy as np
import os

script_dir = os.path.dirname(__file__)
csv_file = os.path.join(script_dir, '../dep/NMR_freq_table.csv')

nuctable=pd.read_csv(csv_file)
gyr_ratio_MHz_T=nuctable["GyrHz"]
name_nuc=nuctable["Name"]

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
        