---
description: Some small applications made using Streamlit and deployed with Render.
---

# NMR Web Apps

## Running the scripts online
Just click on the link beside the script.
It takes a while for the initial run as the service needs to wake up.
As I am using free service of Render, the apps go to sleep when not used for sometime.
Once it starts, everything should be fast.


## Running the scripts offline on a computer
### Cloning the repository / Downloading codes
1. On your computer install git. For windows download it from the [Git for Windows](https://gitforwindows.org/). It is also available on CEA Centre Logiciel.
2. On the git terminal, called Git Bash, clone with the command:
```
git clone https://github.com/subhradip-paul/nmr_webapps.git
````
3. Now you will have a directory called nmr_webapps
4. Navigate to the src folder.
### Preparing Requirements
5. Assuming that you have Python installed in your computer, write:
```
pip install streamlit
```
or run:
````
pip install -r requirements.txt
````
and then run the code using
```
streamlit run <code-name>.py
```
It will start a local server and the script will run from there.


## Summary of scripts
### Dipolar Coupling Calculator [dipole_calculator.py](https://nmrwebapps-dipolecalculator.streamlit.app/)

The script reads the NMR properties of a nucleus from a .csv file and calculates the dipolar coupling based on the input distance in Angstrom. Alternatively it can also calculate the distance if the user inputs the dipolar coupling.

### Temperature Calculator [temp_from_KBr.py](https://temp-from-brt1.onrender.com/)

The script calculates the temperature of a sample, based on the relaxation time of $^79$Br in KBr. More details are in the app.

### Inadequate Efficiency [inadequate_efficiency.py](https://inadequate-efficiency.onrender.com/)

The script calculates the inadequate efficiency based on the relaxation time ($T'_2$) and the J-coupling between the nuclei under observation.

### DNP Sample Prep [dnp_sample_prep.py](https://nmrwebapps-dnpsampleprep.streamlit.app/)

The scripts helps to calculate the *weight* of the biradical for a needed concentration where the user knows the volume needed, or calculate the *volume* of solution for a given weight.

### Maximum sensitivity [max_sensitivity_fromt1.py](https://nmrwebapps-maxsentivityt1.streamlit.app/)

The code gives the recycle delay at which you can ontain the maximum sensitivity. It calculates the build up of the signal as a function of the relaxation time, and then calculates the differential of it. The relaxation modules supported are monoexponential, biexponential, and stretched expoenetial.

### NMR Properties [nmr_properties.py](https://nmr-properties.onrender.com/)
The script gives you the properties of a nuclei and its isotopes, based on the user selection.
