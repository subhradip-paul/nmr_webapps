# Description
Some small applications made with Streamlit and deployed with Render. To run a script you need to install streamlit.<br>*pip install streamlit* <br> and then run the code using <br> *streamlit run <codename>.py*

## NMR Web Apps
<img src=".gitbook/assets/streamlit-mark-color.png" alt="streamlit-logo" data-size="line" width="80" class="center"/>

### Dipolar Coupling Calculator (dipole_calculator.py)

The script reads the NMR properties of a nucleus from a .csv file and calculates the dipolar coupling based on the input distance in Angstrom. Alternatively it can also calculate the distance if the user inputs the dipolar coupling.


### Temperature Calculator (temp_from_KBr.py)

The script calculates the temperature of a sample, based on the relaxation time of $^79$Br in KBr. More details are in the app.

### Inadequate Efficiency (inadequate_efficiency.py)

The script calculates the inadequate efficiency based on the relaxation time ($T'_2$) and the J-coupling between the nuclei under observation.

### DNP Sample Prep (dnp_sample_prep.py)

The scripts helps to calculate the *weight* of the biradical for a needed concentration where the user knows the volume needed, or calculate the *volume* of solution for a given weight.


### Maximum sensitivity (max_sensitivity_fromt1.py)

The code gives the recycle delay at which you can ontain the maximum sensitivity. It calculates the build up of the signal as a function of the relaxation time, and then calculates the differential of it. The relaxation modules supported are monoexponential, biexponential, and stretched expoenetial.

### NMR Properties (nmr_properties.py)
The script gives you the properties of a nuclei and its isotopes, based on the user selection.
