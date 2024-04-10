import streamlit as st

st.set_page_config(
    page_title="Python apps for NMR",
    page_icon=":atom_symbol:",
)

st.title("Python apps for NMR")
st.markdown("Some small apps written by me to help NMR calculations and setting up experiments.\
    The apps are written in [Python](https://www.python.org/) and deployed using [Streamlit](https://streamlit.io/).\
        If you want to help or have suggestions, or if you spot some mistakes, please let me know on the github repository.\
            The github repository is here: https://github.com/subhradip-paul/nmr_webapps.git")
col1, col2 = st.columns(2)
with col1:
    st.image('streamlit-mark-color.png', width=40)
with col2:
    st.image('Python-logo.png', width=40)
st.sidebar.success("Select an app above")