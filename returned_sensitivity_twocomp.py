# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 18:52:52 2022

@author: Paul
"""
#%% import modules
from sympy import symbols, diff, exp, sqrt, lambdify
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

# init_session()

#%% Create the symbolic eqn
t, b, a, Ta, Tb = symbols('t b a Ta Tb')
b = 1 - a
s = a * (1 - exp(-t / Ta)) / sqrt(t) + b * (1 - exp(-t / Tb)) / sqrt(t)
bup = s * sqrt(t)
k = diff(s, t)
st.latex("Sensitivity:" + r'''I_a \times (1 - e^{-t / T_a}) / \sqrt{t} + I_b \times (1 - e^{-t / T_b}) / \sqrt{t}''')

#%% User Inputs
t1=st.number_input("Enter the first buildup time in s: ", value=10.0)
t2 = st.number_input("Enter the second buildup time in s: ", value=100.0)
comp1 = st.number_input("Enter the value of the first component: ", min_value=1.0, value=50.0)
comp2 = st.number_input("Enter the value of the second component: ", min_value=0.1, value=50.0)

comp = comp1 / (comp1 + comp2)


#%% Feed values to the equation
y = k.subs({Ta: t1, Tb: t2, a: comp})
y2 = s.subs({Ta: t1, Tb: t2, a: comp})
y3 = bup.subs({Ta: t1, Tb: t2, a: comp})

lam_x = lambdify(t, y, modules=['numpy'])
lam_x2 = lambdify(t, y2, modules=['numpy'])
lam_x3 = lambdify(t, y3, modules=['numpy'])

x_vals = np.linspace(1, 1024, 100000)
y_vals = lam_x(x_vals)
y2_vals = lam_x2(x_vals)
y3_vals = lam_x3(x_vals)
min_x_val = x_vals[np.argmin(abs(y_vals))]


#%% Plot the equation and the differential
fig, ax = plt.subplots()
ax.plot(x_vals, abs(y_vals), label=r"$\frac{d \left( \frac{I}{\sqrt{t}} \right) }{dt}$")
ax.plot(x_vals, abs(y2_vals), label=r"$\frac{I}{\sqrt{t}}$")
ax.plot(x_vals, y3_vals, label="Intensity(I)")
ax.plot((min_x_val, min_x_val), (0, np.max(y3_vals)), 'k-.')
ax.set_xlabel('Time(s)')
ax.set_ylabel('Intensity (arb. u.)')
ax.legend()
plt.show()

#%% Output of streamlit
st.write('Recycle delay for best sensitivity: ' + str(np.round(x_vals[np.argmin(abs(y_vals))],2)))
st.pyplot(fig)

