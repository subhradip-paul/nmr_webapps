#%% Header Files
import pandas as pd
import numpy as np
import streamlit as st
import os


def main():
    #%% Read the data file
    script_dir = os.path.dirname(__file__)
    csv_file = os.path.join(script_dir, '../dep/biradicals_molarmass.csv')
    df = pd.read_csv(csv_file)
    df = df.sort_values(by=['Name'])
    df2 = df[['Name','MW']]
    st.dataframe(df2, hide_index=True)

    #%% Ask for the biradical and conc.
    rad = st.selectbox('Biradical: ', df.Name, index=1)
    conc = st.number_input('Enter the concentration in mM: ', value = 10.0)
    rad_weight=df["MW"]
    rad_name = df["Name"]
    w_idx = rad_name[rad_name.str.fullmatch(rad)].index

    #%% Ask for the choice of calculation 
    st.markdown('If you have a certain mass of biradical and you want to calculate the volume of solvent \
                you would like to add to make a certain concentration \
                then choose **volume** in the selection below. If you know the volume of the solvent already and \
                would like to know how much biradical you would like to add for a certain concentration then choose **weight** below')
    choice = st.selectbox('Do you want to calculate volume, or weight', ('','volume', 'weight'), index=0)

    #%% main script
    if choice == 'volume':
        w_mg = st.number_input('Weight of biradical in mg: ', value=2.0)
        volume = w_mg*1e6/(conc*rad_weight[w_idx])
        st.write('The volume needed is '+ str(round(volume.iloc[0],2)) + ' $\mu$l')
        solvent=st.selectbox('Choice of solvent: ', ('GDH', 'TCE', 'DMSO'))
        if solvent == 'GDH':
            ratio=0.0
            while ratio != 100.0:
                glycerol=st.number_input('Glycerol percentage: ', value=60.0)
                d2o=st.number_input('D2O percentage: ', value=30.0)
                h2o=st.number_input('H2O percentage: ', value=10.0)
                ratio=glycerol+d2o+h2o
                if ratio != 100.0:
                    st.write('Check percentages')
            st.write('d8-Glycerol needed is ' + str(round(volume.iloc[0]*(glycerol/100)*1.37,2)) + ' mg')
            st.write('$D_2O$ needed is ' + str(round(volume.iloc[0]*(d2o/100),2)) + ' $\mu$l')
            st.write('$H_2O$ needed is ' + str(round(volume.iloc[0]*(h2o/100),2)) + ' $\mu$l')

    if choice == 'weight':
        volume = st.number_input('Volume in $\mu$l: ', value=200.0)
        w_mg=conc*volume*rad_weight[w_idx]/1e6
        w_all=conc*volume*rad_weight/1e6
        st.write('The weight needed is ' + str(round(w_mg.iloc[0],2)) + ' mg')
    
if __name__ == "__main__":
    st.markdown('# DNP Sample Preparation Guide')
    st.markdown('This script is helpful for **DNP** sample preparation \
                The list of biradicals that can be used for calculations are listed below \
                If you would like me to add more biradicals, let me know the name of the biradical and molecular weight\
                If you spot any error or mistakes, email me on subhradip.paul@cea.fr\
                The code comes with absolutely NO WARRANTY, If in doubt, verify the calculations yourself.' )
    st.divider()
    main()        
    