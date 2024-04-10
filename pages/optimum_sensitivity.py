import numpy as np
import matplotlib.pyplot as plt
import streamlit as st


def relax_model(a, ta, tb, beta, xmin, xmax):
  optd1=np.empty([])
  index=np.empty([])
  x=np.linspace(xmin, xmax, 1000000)
  y=1-np.exp(-x/tb)+a*np.exp(-x/tb)-a*np.exp(-(x/ta)**beta)
  sens=y/np.sqrt(x)
  dsens = np.gradient(sens,(xmax-xmin)/1000000)
  index = np.where(np.diff(np.sign(np.squeeze(dsens))) != 0)[0] + 1
  if index.size > 0:
    optd1 = x[index]
  return x, y, sens, optd1


def main():
    option = st.selectbox('Which model?', ('Monoexponential', 'Biexponential', 'Stretched Exponential'))

    if option == 'Monoexponential':
        ta = st.number_input('Enter the time component: ', value = 10.0)
        max_time = st.number_input('Enter the maximum time in the buildup: ', value = 256.0)
        x, snr, sens, optd1 = relax_model(1.0, ta, 1000000.0, 1.0, 0.001, max_time)

    if option == 'Stretched Exponential':
        ta = st.number_input('Enter the time component: ', value = 10.0)
        beta = st.number_input('Enter the $\beta$ component: ', min_value=0.1, max_value=1.00, step=0.01, value = 0.6)
        max_time = st.number_input('Enter the maximum time in the buildup: ', value = 256.0)
        x, snr, sens, optd1 = relax_model(1.0, ta, 1000000.0, beta, 0.001, max_time)

    if option == 'Biexponential':
        ta = st.number_input('Enter the 1st time component: ', value = 1.0)
        tb = st.number_input('Enter the 2nd time component: ', value = 10.0)
        comp1 = st.number_input('Enter \% of 1st component: ', value = 1.0)
        amp1=comp1/100
        max_time = st.number_input('Enter the maximum time in the buildup: ', value = 256.0)
        x, snr, sens, optd1 = relax_model(amp1, ta, tb, 1.0, 0.001, max_time)
    
    fig, ax = plt.subplots()
    ax.plot(x, snr, label='Intensity')
    ax.plot(x, sens, label='Sensitivity')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Intensity')
    ax.legend(loc = 'best')
    ax.set_ylim([0.0, 1.2])
        
    if not np.shape(optd1):
        st.subheader("No suitable solution found, derivative has **NO** zero crossing.")
    else:
        ax.axvline(x=optd1[0], ls = '--', lw = 0.5, c='magenta')
        ax.text(x=optd1[0]+5, y=1.1, s=np.round(optd1[0],2) , fontsize=12)
        st.write('The optimum recycle delay is ' + str(np.round(optd1[0],2)) + ' s')
    
    st.pyplot(fig)
 
if __name__ == "__main__":
    st.title('Optimum recycle delay')
    st.divider()
    st.markdown('The script tries to find the recycle delay at which the sensitivity will be the best.\
                It works for three cases, biexponential buildup, monoexponential buildup, and stretched exponential.\
                It takes the following general equation for buildup:')
    st.latex(r"y = a [1-e^{-(x/T_a)^\beta}] + (1-a) [1-e^{-x/T_b}]")
    st.markdown('In this case sensitivity can be written as:')
    st.latex(r"f(x) = y/\sqrt{x}")
    st.markdown('So the best sensitivity will be when:')
    st.latex(r"\frac{df(x)}{dx} = 0")
    st.write(r"For monoexponential and biexponential, $\beta$ = 1.0")
    st.markdown('For monoexponential and stretched exponential, a = 1.0')
    
    st.divider()
    main()                




