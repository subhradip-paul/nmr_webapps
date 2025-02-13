# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 18:52:52 2022

@author: Paul
"""

#%% Import the modules
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from numpy.polynomial import Polynomial

# #%% Find a value within an array and return the index
# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return array[idx], idx

#%% Find the root

def find_temp_solpol(t1):
    coeff3 = 2.48e09
    coeff2 = 1.42e07
    coeff1 = 5330
    coeff0 = 0.0145 - t1
    p = Polynomial([coeff0,coeff1,coeff2,coeff3])
    k = p.roots()
    temp = np.sqrt(1/k[ k > 0])
    return temp
    



def main():
#%% Create the array for temperature
    
    if 'temp_from_poly' not in st.session_state:
        st.session_state.temp_from_poly = 297.98
    
    temp = np.linspace(10,400,10000)
    y0=0.0145*np.ones(np.size(temp))
    y1=5330*(temp ** -2)
    y2=1.42e7*(temp ** -4)
    y3=2.48e9*(temp ** -6)
    y=y0+y1+y2+y3
    

    #%% Enter the relaxation time and find the temperature
    t1 = st.number_input(":blue[Enter a relaxation time in s:] ", min_value=0.05, max_value=1000.0, step=0.01, value = 0.076333)
    temp_from_poly=find_temp_solpol(t1)[0]
    
    
    

    #%% Create the figure
    fig, ax = plt.subplots()
    ax.plot(temp,np.log10(y))
    ax.axhline(y=np.log10(t1), ls = '--', lw=0.5, c = 'g')
    ax.axvline(x=temp_from_poly, ls = '--', lw=0.5, c = 'r', label=str(round(np.real(temp_from_poly),2))+' K')
    ax.set_xlabel('Relaxation Time (s)')
    ax.set_ylabel('log10 (T) (K)')
    ax.legend(loc='upper right')

    #%% Output cells
    kk = st.session_state.temp_from_poly
    st.metric("**Temperature in K**", np.round(np.real(temp_from_poly),2), delta=np.round(np.real(temp_from_poly-st.session_state.temp_from_poly),2))
    st.session_state.temp_from_poly = temp_from_poly
    
    st.pyplot(fig)


st.title('Temperature from K$$^{79}$$Br')
st.divider()
st.markdown('The temperature is based on the eqn. <br>'  +
            '$$T_1 = 0.0145 + 5330 \\times T^{-2} + 1.42e7 \\times T^{-4} + 2.48e9 \\times T^{-6}$$ <br>' +
            '[Thurber, Tycko. JMR 196, no. 1: 84–87](https://doi.org/10.1016/j.jmr.2008.09.019).'
            , unsafe_allow_html=True)
st.markdown('The equation can be cast as a polynomial substituting $x = T^{-2}$')
st.latex(r"2.49 \times 10^9 x^3 + 1.42 \times 10^7 x^2 + 5330 x + (0.0145 - T) = 0")
st.markdown(r'Solve for x, and $T = \sqrt{1/x}$')
st.divider()
main()
