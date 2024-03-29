import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff, exp, sqrt, lambdify
import streamlit as st

def optimumrecycledelaytwoexp(x, ampl1, tau1, tau2):
    t, b, a, Ta, Tb = symbols ('t b a Ta Tb')
    b = 1-a
    s=a*(1-exp(-t/Ta))/sqrt(t)+b*(1-exp(-t/Tb))/sqrt(t)
    bup=s*sqrt(t)
    k=diff(s,t)
    y_bup=bup.subs({Ta:tau1,Tb:tau2,a:ampl1})
    y=k.subs({Ta:tau1,Tb:tau2,a:ampl1})
    y_sens=s.subs({Ta:tau1,Tb:tau2,a:ampl1})
    lam_x = lambdify(t, y, modules=['numpy'])
    lam_x2 = lambdify(t, y_bup, modules=['numpy'])
    lam_x3 = lambdify(t, y_sens, modules=['numpy'])
    x_vals = np.linspace(0.1, x, 1000000)
    y_vals = lam_x(x_vals)
    y_vals_bup = lam_x2(x_vals)
    y_vals_sens = lam_x3(x_vals)
    optd1diff = x_vals[np.argmin(y_vals)]
    optd1sens = x_vals[np.argmax(y_vals_sens)]
    optd1 = min(optd1diff,optd1sens)                
    return x_vals, y_vals_bup, y_vals_sens, optd1

def optimumrecycledelayoneexp(x, tau1):
    t, Ta = symbols ('t Ta')
    s=(1-exp(-t/Ta))/sqrt(t)
    bup=s*sqrt(t)
    k=diff(s,t)
    y_bup=bup.subs({Ta:tau1})
    y=k.subs({Ta:tau1})
    y_sens=s.subs({Ta:tau1})
    lam_x = lambdify(t, y, modules=['numpy'])
    lam_x2 = lambdify(t, y_bup, modules=['numpy'])
    lam_x3 = lambdify(t, y_sens, modules=['numpy'])
    x_vals = np.linspace(0.1, x, 1000000)
    y_vals = lam_x(x_vals)
    y_vals_bup = lam_x2(x_vals)
    y_vals_sens = lam_x3(x_vals)
    optd1diff = x_vals[np.argmin(y_vals)]
    optd1sens = x_vals[np.argmax(y_vals_sens)]
    optd1 = min(optd1diff,optd1sens)   
    return x_vals, y_vals_bup, y_vals_sens, optd1

def optimumrecycledelaystrexp(x, tau1, beta):
    t, b, Ta = symbols ('t b Ta')
    s=(1-exp(-(t/Ta)**b))/sqrt(t)
    bup=s*sqrt(t)
    k=diff(s,t)
    y_bup=bup.subs({Ta:tau1,b:beta})
    y=k.subs({Ta:tau1,b:beta})
    y_sens=s.subs({Ta:tau1,b:beta})
    lam_x = lambdify(t, y, modules=['numpy'])
    lam_x2 = lambdify(t, y_bup, modules=['numpy'])
    lam_x3 = lambdify(t, y_sens, modules=['numpy'])
    x_vals = np.linspace(0.1, x, 1000000)
    y_vals = lam_x(x_vals)
    y_vals_bup = lam_x2(x_vals)
    y_vals_sens = lam_x3(x_vals)
    optd1diff = x_vals[np.argmin(y_vals)]
    optd1sens = x_vals[np.argmax(y_vals_sens)]
    optd1 = min(optd1diff,optd1sens)                
    return x_vals, y_vals_bup, y_vals_sens, optd1

option = st.selectbox('Which model?', ('Monoexponential', 'Biexponential', 'Stretched Exponential'))


if option == 'Monoexponential':
    tb = st.number_input('Enter the time component: ', value = 10.0)
    max_time = st.number_input('Enter the maximum time in the buildup: ', value = 256.0)
    x, snr, sens, optd1 = optimumrecycledelayoneexp(max_time, tb)

if option == 'Stretched Exponential':
    tb = st.number_input('Enter the time component: ', value = 10.0)
    beta = st.number_input('Enter the $\beta$ component: ', min_value=0.55, max_value=1.00, step=0.01, value = 0.6)
    max_time = st.number_input('Enter the maximum time in the buildup: ', value = 256.0)
    x, snr, sens, optd1 = optimumrecycledelaystrexp(max_time, tb, beta)

if option == 'Biexponential':
    ta = st.number_input('Enter the 1st time component: ', value = 1.0)
    tb = st.number_input('Enter the 2nd time component: ', value = 10.0)
    comp1 = st.number_input('Enter the Amplitude 1: ', value = 1.0)
    comp2 = st.number_input('Enter the Amplitude 1: ', value = 99.0)
    amp1=comp1/(comp1+comp2)
    max_time = st.number_input('Enter the maximum time in the buildup: ', value = 256.0)
    x, snr, sens, optd1 = optimumrecycledelaytwoexp(max_time, amp1, ta, tb)
    

fig, ax = plt.subplots()
ax.plot(x, snr, label='SNR')
ax.plot(x, sens, label='Sensitivity')
ax.plot((optd1,optd1),(0.0,1.0),'k--')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Intensity')
ax.legend(loc = 'best')
ax.set_ylim([0.0, 1.2])

st.write('The optimum recycle delay is ' + str(np.round(optd1,2)) + ' s')
st.pyplot(fig)
