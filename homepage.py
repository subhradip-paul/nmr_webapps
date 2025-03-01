import streamlit as st


home = st.Page('_pages/home.py', default=True, title='Home Page', icon='🔬')
nuclear_properties = st.Page("codes/NMR_Nuclei_Parameters.py", title="Nuclear properties", icon="⚛")
optimum_recycle_delay = st.Page("codes/Optimum_Recycle_Delay.py", title="Optimum Recycle Delay", icon="📈")
sample_temp_kbr = st.Page("codes/Sample_Temp_from_KBr_T1.py", title ="Sample Temperature from KBr T1", icon="🌡️")
dnp_sample_prep = st.Page("codes/DNP_Sample_Preparation.py", title ="DNP Sample Preparation", icon="🍲")
dipole_calculator = st.Page("codes/Dipole_Calculator.py", title="Dipolar Coupling Calculator", icon="↔")
inadequate_efficiency_calc = st.Page("codes/INADEQUATE_Efficiency.py", title="Setting up INADEQUATE", icon ="⚖️")
power_calculator = st.Page("codes/Power_Calculator.py", title ="Power and pulse lengths in NMR", icon="⚡")
biradical_properties = st.Page("codes/Biradical_Visualiser.py", title = "Structure of DNP Radicals", icon="🧬")


pg = st.navigation(
{
        "Home": [home],
        "Nuclei in NMR": [nuclear_properties],
        "Setting up Experiments": [optimum_recycle_delay, sample_temp_kbr, dipole_calculator, inadequate_efficiency_calc, power_calculator],
        "DNP Related": [dnp_sample_prep, biradical_properties],
    }
)

pg.run()