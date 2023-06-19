#%% Header files
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import plotly.graph_objects as go


#%% Plotting function for inadequate
def plotinad(t2prime, jcc):
    x=np.linspace(0,50,100000)
    y=(np.sin(2*np.pi*jcc*x*1e-3)**2)*np.exp(-4*x/t2prime)
    fig = go.Figure(data=go.Scatter(x=x, y=y))
    fig.update_layout(title = "Inadequate DQ Efficiency", 
                      xaxis_title=r"$ \tau_0 $ (ms)", yaxis_title='<b> DQ Efficiency </b>',
                      font_size = 28,
                      font_family = 'Lucida Fax')
    return fig

#%% User Input
st.markdown(r"This is a plot of DQ Efficiency of INADEQUATE, based \
    on the paper Cadars et al. *J. Phys. Chem. B.* 110 (2006) 16982."    )

t2prime_input=st.number_input("$T'_2$ (ms)", min_value=0.1, max_value=200.0, value=10.0, step=0.5)
Jcc_input=st.number_input("$J_{II}$ (Hz)", min_value=1.0, max_value=500.0, value=54.0, step=0.5)

st.markdown(r"**The total mixing time is 4 $\tau_0$**"    )

fig_inad=plotinad(t2prime_input, Jcc_input)


st.plotly_chart(fig_inad)