# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 18:52:52 2022

@author: Paul
"""

#%% Import the modules
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

#%% Find a value within an array and return the index
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def main():
#%% Create the array for temperature
    st.title('Temperature from K$$^{79}$$Br')
    st.subheader('Created by Subhradip Paul')
    st.markdown('The temperature is based on the eqn. <br>'  +
                '$$T_1 = 0.0145 + 5330 \\times T^{-2} + 1.42e7 \\times T^{-4} + 2.48e9 \\times T^{-6}$$ <br>' +
                '[Thurber, Tycko. JMR 196, no. 1: 84–87](https://doi.org/10.1016/j.jmr.2008.09.019).' 
                , unsafe_allow_html=True)
    temp = np.linspace(20,400,4000)
    y0=0.0145*np.ones(np.size(temp));
    y1=5330*(temp ** -2);
    y2=1.42e7*(temp ** -4);
    y3=2.48e9*(temp ** -6);
    y=y0+y1+y2+y3;
    

    #%% Enter the relaxation time and find the temperature
    t1 = st.number_input("Enter a relaxation time in s: ", value = 15)
    k, l=find_nearest(y, t1)


    #%% Create the figure
    fig, ax = plt.subplots()
    ax.plot(temp,y)
    ax.plot((temp[1], temp[np.size(temp)-1]), (t1, t1), 'k--')
    ax.set_xlabel('Relaxation Time (s)')
    ax.set_ylabel('Temperature (K)')

    #%% Output cells
    st.write('The temperature is '+ str(np.round(temp[l],2)) + ' K')
    st.pyplot(fig)

if __name__ == "__main__":
    main()