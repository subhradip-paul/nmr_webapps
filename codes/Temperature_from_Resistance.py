#!/usr/bin/env python3
import numpy as np
import streamlit as st
from numpy.polynomial import Polynomial





class ResistanceToTemperatureConverter:
    def __init__(self):        pass
    
    def  CX95572(self, resistance_ohm):
        if resistance_ohm > 1187:
            coefficients = [-8.66315e-2, 8.51080e-2, -1.97582e-2, 1.38816e-3]
        else:
            coefficients = [-1.15282e-1, 6.63234e-2, -1.29163e-2, 8.77272e-4]
        inverse_temperature = sum(
            coefficient * np.log(resistance_ohm) ** power
            for power, coefficient in enumerate(coefficients)
        )
        return 1 / inverse_temperature
    
    def PTN(self, resistance_ohm, n=1):
        r0 = 100.0 * n
        A = 3.9083e-3
        B = -5.775e-7
        C = -4.183e-12

        def resistance_at(T_k):
            t = T_k - 273.15  # °C
            if t >= 0:
                return r0 * (1 + A*t + B*t*t)
            return r0 * (1 + A*t + B*t*t + C*(t - 100)*t**3)

        # Nominal IEC 60751 operating interval: 73.15–298 K
        lo, hi = 73.15, 298.0
        r_lo, r_hi = resistance_at(lo), resistance_at(hi)

        if not r_lo <= resistance_ohm <= r_hi:
            raise ValueError(
                f"Resistance {resistance_ohm:.4f} Ω is outside the "
                f"{lo}–{hi} K nominal range ({r_lo:.4f}–{r_hi:.4f} Ω)."
            )

        # Bisection: converge to the temperature whose calculated R matches input R
        for _ in range(50):
            mid = (lo + hi) / 2
            if resistance_at(mid) < resistance_ohm:
                lo = mid
            else:
                hi = mid

        return (lo + hi) / 2


def main():
    st.title("Cernox Temperature Converter")

    resistance = st.number_input(
        "Enter resistance (Ω):",
        min_value=0.001,  # avoids log(0)
        value=166.0,
    )
    
    sensor_type = st.selectbox(
        "Select sensor type:",
        ["CERNOX X95572 (Sample)", "PTN (Sample)"],
    )

    if sensor_type is "CERNOX X95572 (Sample)":
        temperature = ResistanceToTemperatureConverter().CX95572(resistance)
    elif sensor_type is "PTN (Sample)":
        n = st.number_input(
            "Enter n value for PTN sensor (default is 1):",
            min_value=1,
            max_value=1000,
            value=100,
        )
        temperature = ResistanceToTemperatureConverter().PTN(resistance, n)
    st.write(f"Estimated temperature: {temperature:.2f} K")


if __name__ == "__main__":
    main()