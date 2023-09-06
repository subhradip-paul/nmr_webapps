---
description: Some small applications made with flask or Streamlit.
---

# NMR Web Apps

Some small applications made with Streamlit and deployed with Render. To run a script you need to install streamlit.\
_pip install streamlit_\
and then run the code using\
_streamlit run .py_

## NMR Web Apps

<img src=".gitbook/assets/streamlit-mark-color.png" alt="streamlit-logo" data-size="line">

### Dipolar Coupling Calculator (dipole\_calculator.py)

The script reads the NMR properties of a nucleus from a .csv file and calculates the dipolar coupling based on the input distance in Angstrom. Alternatively it can also calculate the distance if the user inputs the dipolar coupling.

### Temperature Calculator (temp\_from\_KBr.py)

The script calculates the temperature of a sample, based on the relaxation time of $^79$Br in KBr. More details are in the app.



D

### Inadequate Efficiency (inadequate\_efficiency.py)

The script calculates the inadequate efficiency based on the relaxation time ($T'\_2$) and the J-coupling between the nuclei under observation.

### DNP Sample Prep (dnp\_sample\_prep.py)

The scripts helps to calculate the _weight_ of the biradical for a needed concentration where the user knows the volume needed, or calculate the _volume_ of solution for a given weight.

### Maximum sensitivity (max\_sensitivity\_fromt1.py)

The code gives the recycle delay at which you can ontain the maximum sensitivity. It calculates the build up of the signal as a function of the relaxation time, and then calculates the differential of it. The relaxation modules supported are monoexponential, biexponential, and stretched expoenetial.

### NMR Properties (nmr\_properties.py)

The script gives you the properties of a nuclei and its isotopes, based on the user selection.
