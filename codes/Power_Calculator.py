import numpy as np
import streamlit as st

def kHz2watts(pow1, pow2, watt1):
    watt2 = watt1 * (pow2/pow1)**2
    return watt2

def plen2watts(plen1, plen2, watt1):
    watt2 = kHz2watts(plen2, plen1, watt1)
    return watt2

def plen2db(plen1, plen2, db1):
    deltadb = 20 * np.log10(plen1/plen2)
    db2 = deltadb + db1
    return db2

def db2plen(db1, db2, plen1):
    deltadb = db1-db2
    plen2 = plen1 * (10 ** (deltadb/20))
    return plen2

def main():
    choice = st.selectbox('What do you want to calculate?', (' ', 'Power in Watts', 'Pulse length for a given Watt', 'Power in dB', 'Pulse length for a given dB'))

    if choice == 'Power in Watts':
        pow_ini=st.number_input('Initial Power in kHz :', value=100)
        pow_final=st.number_input('Final Power in kHz :', value=50)
        pow_watts=st.number_input('Initial Power in Watts:', value=70)
        pow_final=kHz2watts(pow_ini, pow_final, pow_watts)
        st.metric('Power needed in Watts', pow_final, (pow_final-pow_ini))
    elif choice == 'Pulse length for a given Watt':
        plen_ini=st.number_input(r'Initial Pulse Length in $\mu$s :', value=2.5)
        plen_final=st.number_input(r'Final Pulse Length in $\mu$s :', value=5.)
        pow_watts=st.number_input('Initial Power in Watts:', value=70)
        pow_final=plen2watts(plen_ini, plen_final, pow_watts)
        st.metric('Power needed in Watts', pow_final, (plen_ini-plen_final))
    elif choice == 'Power in dB':
        plen_ini=st.number_input(r'Initial Pulse Length in $\mu$s :', value=2.5)
        plen_final=st.number_input(r'Final Pulse Length in $\mu$s :', value=5.)
        pow_db=st.number_input('Initial Power in dB:', value=70)
        pow_final=plen2db(plen_ini, plen_final, pow_db)
        powchange=np.round((pow_final-pow_db), 2)
        st.metric('Power needed in dB', np.round(pow_final,2), powchange)
    elif choice == 'Pulse length for a given dB':
        plen_ini=st.number_input(r'Initial Pulse Length in $\mu$s :', value=2.5)
        pow_db1=st.number_input(r'Initial dB :', value=70)
        pow_db2=st.number_input('Final dB:', value=64)
        plen_final=db2plen(pow_db1, pow_db2, plen_ini)
        st.metric(r"Final pulse in $\mu$s", np.round(plen_final,2), np.round(plen_final-plen_ini,2))
        
st.title('Power Calculator')
st.divider()
main()