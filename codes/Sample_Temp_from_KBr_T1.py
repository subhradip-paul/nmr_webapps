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
import pandas as pd



#%% Find the root

def find_temp_solpol(t1, c1, c2, c3, c4):
    coeff3 = c4
    coeff2 = c3
    coeff1 = c2
    coeff0 = c1 - t1
    p = Polynomial([coeff0,coeff1,coeff2,coeff3])
    k = p.roots()
    temp = np.sqrt(1/k[ k > 0])
    return temp
    



def main():
#%% Create the array for temperature
    sample = st.selectbox("Choose sample: ", options = ["KBr", "CsBr", "CsI", "RbBr"])

    if 'temp_from_poly' not in st.session_state:
        st.session_state.temp_from_poly = 297.98


    coeff_1 = 0.0145
    coeff_2 = 5.33e3
    coeff_3 = 14.2e6
    coeff_4 = 2.48e9

    if sample == "KBr":
        coeff_1 = 0.0145
        coeff_2 = 5.33e3
        coeff_3 = 14.2e6
        coeff_4 = 2.48e9
    elif sample == "RbBr":
        coeff_1 = 5.9e-3
        coeff_2 = 4.74e3
        coeff_3 = 5.9e6
        coeff_4 = 2.08e9
    elif sample == "CsBr":
        coeff_1 = 5.5e-3
        coeff_2 = 6.47e3
        coeff_3 = 3.94e6
        coeff_4 = 0.52e9
    elif sample == "CsI":
        coeff_1 = -1.6e-3
        coeff_2 = 1.52e3
        coeff_3 = 0.387e6
        coeff_4 = 0.121e9





    temp = np.linspace(10,400,10000)
    y0=coeff_1*np.ones(np.size(temp))
    y1=coeff_2*(temp ** -2)
    y2=coeff_3*(temp ** -4)
    y3=coeff_4*(temp ** -6)
    y=np.real(y0+y1+y2+y3)
    

    #%% Enter the relaxation time and find the temperature
    t1 = st.number_input(":blue[Enter a relaxation time in s:] ", min_value=0.05, max_value=1000.0, step=0.01, value = 0.076333)
    temp_from_poly=find_temp_solpol(t1, coeff_1, coeff_2, coeff_3, coeff_4)[0]
    
    
    

    #%% Create the figure
    # fig, ax = plt.subplots()
    # ax.plot(temp,np.log10(y))
    # ax.axhline(y=np.log10(t1), ls = '--', lw=0.5, c = 'g')
    # ax.axvline(x=temp_from_poly, ls = '--', lw=0.5, c = 'r', label=str(round(np.real(temp_from_poly),2))+' K')
    # ax.set_xlabel('Relaxation Time (s)')
    # ax.set_ylabel('log10 (T) (K)')
    # ax.legend(loc='upper right')
    import plotly.graph_objects as go


    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=temp,
        y=np.log10(np.real(y)),
        mode='lines',
        name='log10(y)'
    ))

    fig.add_hline(
        y=np.log10(t1),
        line_dash='dash',
        line_color='green',
        annotation_text='T1'
    )

    fig.add_vline(
        x=np.real(temp_from_poly),
        line_dash='dash',
        line_color='red',
        annotation_text=f'{np.real(temp_from_poly):.2f} K'
    )

    fig.update_layout(
        xaxis_title="Relaxation Time (s)",
        yaxis_title=r"log{_10} (T) (K)",
        legend=dict(x=0.75, y=0.95)
    )

    # ===== Styling =====
    fig.update_layout(
        xaxis=dict(
            title="Relaxation Time (s)",
            showline=True,
            linewidth=0.75,
            linecolor="black",
            mirror=True,
            showgrid=True,
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False,
        ),
        yaxis=dict(
            title=r"log10 (T) (K)",
            showline=True,
            linewidth=0.75,
            linecolor="black",
            mirror=True,
            showgrid=True,
            gridcolor="rgba(0,0,0,0.08)",
            zeroline=False
        ),
        # background shade
        plot_bgcolor="rgba(240,240,240,0.6)",   # light grey
        paper_bgcolor="white",
        legend=dict(x=0.75, y=0.95),
        margin=dict(l=60, r=20, t=40, b=60)
    )

    st.plotly_chart(fig, width='stretch')

    #%% Output cells
    kk = st.session_state.temp_from_poly
    st.metric("**Temperature in K**", np.round(np.real(temp_from_poly),2), delta=np.round(np.real(temp_from_poly-st.session_state.temp_from_poly),2))
    st.session_state.temp_from_poly = temp_from_poly
    
    # st.pyplot(fig)


st.title(r"Sample temperature from T$$_{1n}$$")
st.divider()
st.markdown('The temperature is based on the eqn. <br>'  +
            '$$T_1 = c_0 + c_2 \\times T^{-2} + c_4 \\times T^{-4} + c_6 \\times T^{-6}$$ <br>' +
            '[Thurber, Tycko. JMR 196 (2009) 84](https://doi.org/10.1016/j.jmr.2008.09.019). <br>' +
            '[Sarkar et al., JMR 212 (2011) 460](https://doi.org/10.1016/j.jmr.2011.08.021)'
            , unsafe_allow_html=True)
st.markdown("The coefficient for different samples are:" )
coeff_data = {'Sample': ['KBr','RbBr','CsBr','CsI'],
              'Nucleus': [r"79 Br", '79 Br', '79 Br', '127 I'],
              'c_0': [14.5e-3, 5.9e-3, 5.5e-3, -1.6e-3],
              'c_2':[5.33e3, 4.74e3, 6.47e3, 1.52e3],
              'c_4': [14.2e6, 5.9e6, 3.94e6, 0.387e6],
              'c_6':[2.48e9, 2.08e9, 0.52e9, 0.121e9]}
data_df = pd.DataFrame(coeff_data)
st.dataframe(data_df, hide_index=True,
             column_config={'c_0': st.column_config.NumberColumn(r"c_0", format="%.2e"),
                                                      'c_2': st.column_config.NumberColumn(r"c_2", format="%.2e"),
                                                      'c_4': st.column_config.NumberColumn(r"c_4", format="%.2e"),
                                                      'c_6': st.column_config.NumberColumn(r"c_6", format="%.2e")})
st.markdown('The equation can be cast as a polynomial substituting $x = T^{-2}$')
st.latex(r"c_6 x^4 + c_4 x^2 + c_2 x + (c_0 - T) = 0")
st.markdown(r'Solve for x, and $T = \sqrt{1/x}$')
st.divider()
main()
