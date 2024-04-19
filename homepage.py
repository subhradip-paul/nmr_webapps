import streamlit as st
from st_pages import Page, Section, show_pages, add_indentation

st.image('dnp_grenoble_logo.png')

st.markdown("Some small apps written by me to help NMR calculations and setting up experiments.\
    The apps are written in [Python](https://www.python.org/) and deployed using [Streamlit](https://streamlit.io/).\
        If you want to help or have suggestions, or if you spot some mistakes, please let me know on the github repository.\
            The github repository is here: https://github.com/subhradip-paul/nmr_webapps.git")
col1, col2 = st.columns(2)
with col1:
    st.image('streamlit-mark-color.png', width=40)
with col2:
    st.image('Python-logo.png', width=40)

add_indentation()
show_pages(
    [   
        Page("homepage.py", "Home", ":house:"),
        Section("Nuclear Properties", icon=":atom_symbol:"),
        Page("pages/nmr_properties.py", "NMR Properties"),
        Page("pages/dipole_calculator.py", "Distance to Dipole and vice versa"),
        Section("Laboratory Stuff", icon="🥼"),
        Page("pages/dnp_sample_prep.py", "DNP Sample Prep"),
        Page("pages/power_calculator.py", "Power Calculator"),
        Section("Experiments", icon="🧨"),
        Page("pages/inadequate_efficiency.py", "INADEQUATE Efficiency"),
        Page("pages/optimum_sensitivity.py", "Optimum Sensitivity based on T1"),
        Page("pages/temp_from_KBr.py", "Temperature from KBr")
    ]
)