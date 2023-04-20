#%% Header Files
import pandas as pd
import numpy as np
import streamlit as st

#%% Read the data file
df = pd.read_csv('biradicals_molarmass.csv')

#%% Ask for the biradical and conc.
rad = st.selectbox('Biradical: ', df.Name)
conc = st.number_input('Enter the concentration in mM: ', value = 10.0)
rad_weight=df["MW"]
rad_name = df["Name"]
w_idx = rad_name[rad_name.str.fullmatch(rad)].index

st.markdown(Please)
choice = st.selectbox('Do you ')