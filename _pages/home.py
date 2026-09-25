import streamlit as st



st.set_page_config(page_title="Home", page_icon="🔬")

st.image('dnp_grenoble_logo.png')
col1, col2 = st.columns(2, gap="small")
with col1:
    st.image('streamlit-mark-color.png', width=30)
with col2:
    st.image('Python-logo.png', width=30)
st.markdown("""
---
### NMR-related calculators and tools for the lab, built with Streamlit as a single multi-page app.
---

# NMR Web Apps

## Running the app online
The apps are bundled into one multi-page Streamlit site — use the sidebar to jump between tools.

**Live app:** [Online App](https://nmr-webapps-28ab810725dd.herokuapp.com)

It can take a moment to load on the first visit if the service has been asleep; after that it's fast.

## Running it offline on your computer
### 1. Get the code
Install [git](https://gitforwindows.org/) (also available on CEA Centre Logiciel), then:
```
git clone https://github.com/subhradip-paul/nmr_webapps.git
cd nmr_webapps
```

### 2. Install dependencies
This project's dependencies are defined in `pyproject.toml` (locked with `uv.lock`). Either:
```
uv sync
```
or, with plain pip:
```
pip install .
```

### 3. Run it
```
streamlit run homepage.py
```
This starts a local server and opens the same multi-page app you'd see online.

## Summary of tools
The sidebar groups tools the same way as below.

### Nuclei and Interactions in NMR
- **Nuclear properties** — [`codes/NMR_Nuclei_Parameters.py`](codes/NMR_Nuclei_Parameters.py): look up NMR properties of a nucleus and its isotopes.
- **Dipolar Coupling Calculator** — [`codes/Dipole_Calculator.py`](codes/Dipole_Calculator.py): calculates the dipolar coupling from an internuclear distance (in Å), or the distance from a given coupling.
- **Chemical Shift Predictor** — [`codes/Chemical_Shift_Prediction.py`](codes/Chemical_Shift_Prediction.py): predicts/looks up chemical shifts, pulling comparable compounds from PubChem.

### Setting up Experiments
- **Optimum Recycle Delay** — [`codes/Optimum_Recycle_Delay.py`](codes/Optimum_Recycle_Delay.py): finds the recycle delay for maximum sensitivity from a relaxation-time build-up curve (mono-exponential, bi-exponential, and stretched-exponential models supported).
- **Sample T from T1n** — [`codes/Sample_Temp_from_KBr_T1.py`](codes/Sample_Temp_from_KBr_T1.py): calculates sample temperature from the ⁷⁹Br T₁ relaxation time in KBr.
- **Setting up INADEQUATE** — [`codes/INADEQUATE_Efficiency.py`](codes/INADEQUATE_Efficiency.py): calculates refocussed INADEQUATE efficiency from T₂′ and the J-coupling.
- **DQ SQ Spectrum Generator** — [`codes/DQ_SQ_predictor.py`](codes/DQ_SQ_predictor.py): generates predicted DQ-SQ spectra.
- **Cernox Temperature Converter** — [`codes/Temperature_from_Resistance.py`](codes/Temperature_from_Resistance.py): converts Cernox sensor resistance to temperature.

### DNP Related
- **DNP Sample Preparation** — [`codes/DNP_Sample_Preparation.py`](codes/DNP_Sample_Preparation.py): calculates the biradical weight needed for a target concentration/volume, or vice versa.
- **Structure of DNP Radicals** — [`codes/Biradical_Visualiser.py`](codes/Biradical_Visualiser.py): looks up and visualizes biradical structures from SMILES data.

### Lab Related
- **Lab Members** — [`codes/Lab_Members.py`](codes/Lab_Members.py): displays lab member info, sourced from a published Google Sheet.
- **Citation Map Generator** — [`codes/Citation_Map_Generator.py`](codes/Citation_Map_Generator.py): builds a citation network map from OpenAlex data.

### Miscellaneous
- **Molecular Weight Calculator** — [`codes/Molecular_Weight.py`](codes/Molecular_Weight.py): calculates molecular weight from a chemical formula.

""")
