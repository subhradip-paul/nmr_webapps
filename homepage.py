import streamlit as st


home = st.Page('_pages/home.py', default=True, title='Home Page', icon='🔬')
nuclear_properties = st.Page("codes/NMR_Nuclei_Parameters.py", title="Nuclear properties", icon="⚛")
optimum_recycle_delay = st.Page("codes/Optimum_Recycle_Delay.py", title="Optimum Recycle Delay", icon="📈")
sample_temp_kbr = st.Page("codes/Sample_Temp_from_KBr_T1.py", title =r"Sample T from T1n", icon="🌡️")
dnp_sample_prep = st.Page("codes/DNP_Sample_Preparation.py", title ="DNP Sample Preparation", icon="🍲")
dipole_calculator = st.Page("codes/Dipole_Calculator.py", title="Dipolar Coupling Calculator", icon="↔")
nmr_shift_predictor = st.Page("codes/Chemical_Shift_Prediction.py", title = "Chemical Shift Predictor", icon = "🤺")
inadequate_efficiency_calc = st.Page("codes/INADEQUATE_Efficiency.py", title="Setting up INADEQUATE", icon ="⚖️")
biradical_properties = st.Page("codes/Biradical_Visualiser.py", title = "Structure of DNP Radicals", icon="🧬")
dq_sq_prediction = st.Page("codes/DQ_SQ_predictor.py", title = "DQ SQ Spectrum Generator", icon="🏋️‍♂️")


pg = st.navigation(
{
        "Home": [home],
        "Nuclei and Interactions in NMR": [nuclear_properties, dipole_calculator, nmr_shift_predictor],
        "Setting up Experiments": [optimum_recycle_delay, sample_temp_kbr, inadequate_efficiency_calc, dq_sq_prediction],
        "DNP Related": [dnp_sample_prep, biradical_properties],
    }
)

pg.run()