#%% Header files
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import streamlit.components.v1 as components
import matplotlib.pyplot as plt, mpld3
from mpld3 import plugins


#%% Plotting function for inadequate
def plotinad(t2prime, jcc):
    x=np.linspace(0,100,100000)
    y=(np.sin(2*np.pi*jcc*x*1e-3)**2)*np.exp(-4*x/t2prime)
    fig, ax = plt.subplots()
    ax = plt.plot(x, y)
    plt.xlabel('Time (ms)')
    plt.ylabel('DQ Efficiency')
    return fig

#%% User Input
st.markdown("This is a plot of DQ Efficiency of INADEQUATE, based \
    on the paper Cadars et al. *J. Phys. Chem. B.* 110 (2006) 16982")

t2prime_input=st.slider("$T'_2$ (ms)", min_value=0.1, max_value=200.0, value=10.0, step=0.5)
Jcc_input=st.slider("$J_{cc}$ (Hz)", min_value=1.0, max_value=100.0, value=54.0, step=0.5)

fig_inad=plotinad(t2prime_input, Jcc_input)

# #st.pyplot(fig_inad)
# fig_html = mpld3.fig_to_html(fig_inad)
# components.html(fig_html, height=600)

# %% Making Interaction

# Define some CSS to control our custom labels
css = '''
table
{
  border-collapse: collapse;
}
th
{
  color: #ffffff;
  background-color: #000000;
}
td
{
  background-color: #cccccc;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 1px solid black;
  text-align: right;
}
'''

for axes in fig_inad.axes:
  for line in axes.get_lines():
    xy_data = line.get_xydata()
    labels = []
    for x,y in xy_data:
      html_label = f'<table border="1" class="dataframe"> <thead> <tr style="text-align: right;"> </thead> <tbody> <tr> <th>x</th> <td>{x}</td> </tr> <tr> <th>y</th> <td>{y}</td> </tr> </tbody> </table>'
      labels.append(html_label)
    tooltip = plugins.PointHTMLTooltip(line, labels, css=css)
    plugins.connect(fig_inad, tooltip)

fig_html = mpld3.fig_to_html(fig_inad)
components.html(fig_html, height=850, width=850)