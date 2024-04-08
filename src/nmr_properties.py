import streamlit as st
import pandas as pd
import numpy as np
import os
import re

script_dir = os.path.dirname(__file__)
datafile = os.path.join(script_dir, '../dep/IsotopeProperties.csv')

df=pd.read_csv(datafile)

st.title("NMR Properties :atom_symbol:")


nuc=st.text_input('Please enter a nucleus name in the format like 1H or H: ', value = "1H")
mag_field=st.number_input('Enter a magnetic field in T: ', value = 9.4)

user_input=re.split('(\d+|[A-Za-z]+)',nuc)
if (len(user_input) == 3):
    nuc_choice=df.loc[df['Nucleus'] == user_input[1]]
    st.markdown(nuc_choice.iloc[:,:10].to_markdown(index=False))
    st.caption("Properties of the Nucleus, Database from [ssNake](https://github.com/smeerten/ssnake).")
    larmor_freq=nuc_choice["Gamma (MHz/T)"].to_numpy().astype(float)
    st.write(pd.DataFrame({
    'Mass Num': nuc_choice["Mass"],
    'Nucleus' : nuc_choice["Nucleus"],
    'Larmor Freq (MHz)': larmor_freq*mag_field,
    }).to_markdown(index=False))
    st.caption("Larmor Frequency :tornado:")
else:
    nuc_choice = df.loc[df['Nucleus'] == user_input[3]] 
    mass_choice = nuc_choice.loc[nuc_choice['Mass'] == user_input[1]]
    st.markdown(mass_choice.iloc[:,:10].to_markdown(index=False))
    st.caption("Properties of the Nucleus, Database from [ssNake](https://github.com/smeerten/ssnake).")
    larmor_freq=mass_choice["Gamma (MHz/T)"].to_numpy().astype(float)
    st.write(pd.DataFrame({
    'Mass Num': mass_choice["Mass"],
    'Nucleus' : mass_choice["Nucleus"],
    'Larmor Freq (MHz)': larmor_freq*mag_field,
    }).to_markdown(index=False))
    st.caption("Larmor Frequency :tornado:")

    
