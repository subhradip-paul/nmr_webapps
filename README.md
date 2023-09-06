---
description: Some small applications made with Streamlit and deployed with Render.
---

# NMR Web Apps

### Dipolar Coupling Calculator (dipole_calculator.py)

The script reads the NMR properties of a nucleus from a .csv file and calculates the dipolar coupling based on the input distance in Angstrom. Alternatively it can also calculate the distance if the user inputs the dipolar coupling.

$$d(r_{12}) = - \frac{\mu_0 \gamma_1 \gamma_2\hbar}{4\pi r_{12}^3}$$

### Temperature Calculator (temp_from_KBr.py)

The script calculates the temperature of a sample, based on the relaxation time of $^79$Br in KBr. More details are in the app.

### Inadequate Efficiency (inadequate_efficiency.py)

The script calculates the inadequate efficiency based on the relaxation time ($$T'_2$$) and the J-coupling between the nuclei under observation.


<img src=".gitbook/assets/streamlit-mark-color.png" alt="" data-size="line">
