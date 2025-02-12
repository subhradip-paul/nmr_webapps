#%% Header files
import numpy as np
import streamlit as st
import plotly.graph_objects as go


#%% Plotting function for inadequate
def plotinad(t2prime, jcc):
    x=np.linspace(0,50,100000)
    inad_eff=(np.sin(2*np.pi*jcc*x*1e-3)**2)*np.exp(-4*x/t2prime)
    fig = go.Figure(data=go.Scatter(x=x, y=inad_eff))
    fig.add_vline(x=x[np.argmax(inad_eff)], line_width=0.5, line_dash="dash", line_color="green")
    fig.update_layout(title = "Inadequate DQ Efficiency", 
                      xaxis_title= "T<sub>0</sub> (ms)", yaxis_title='<b> DQ Efficiency </b>',
                      font_size = 16,
                      font_family = 'Lucida Fax'
                      )
    return fig, inad_eff

def main():
    t2prime_input=st.number_input("$T'_2$ (ms)", min_value=0.1, max_value=200.0, value=10.0, step=0.5)
    Jcc_input=st.number_input("$J_{II}$ (Hz)", min_value=1.0, max_value=500.0, value=54.0, step=0.5)
    x=np.linspace(0,50,100000)
    fig_inad, inad_eff=plotinad(t2prime_input, Jcc_input)
    st.plotly_chart(fig_inad)
    st.write('The max efficiency is: ' + str(round(np.max(inad_eff),2)) + r', and $ \tau_0 $ = ' + str(round(x[np.argmax(inad_eff)],2))+ ' ms')
    st.markdown(r"**The total mixing time is 4 $\tau_0$**"    )
    
if __name__ == "__main__":
    st.title("Refocussed INADEQUATE Efficiency")
    st.divider()
    st.markdown(r"This is a plot of DQ Efficiency of INADEQUATE, based \
    on the papers \
        1. Nakai, T., and C. A. Mcdowell. *J. Magn. Reson.* 104 (1993) 146\
        2. Cadars et al. *J. Phys. Chem. B.* 110 (2006) 16982."    )
    main()