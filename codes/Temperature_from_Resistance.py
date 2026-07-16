#!/usr/bin/env python3
import numpy as np
import streamlit as st


class CernoxCalibration:
    def __init__(self, coefficients):
        self.coefficients = np.asarray(coefficients, dtype=float)

    def temperature_from_resistance(self, resistance_ohm):
        inverse_temperature = sum(
            coefficient * np.log(resistance_ohm) ** power
            for power, coefficient in enumerate(self.coefficients)
        )
        return 1 / inverse_temperature


def main():
    st.title("Cernox Temperature Converter")

    resistance = st.number_input(
        "Enter resistance (Ω):",
        min_value=0.001,  # avoids log(0)
        value=166.0,
    )

    if resistance > 1187:
        coefficients = [-8.66315e-2, 8.51080e-2, -1.97582e-2, 1.38816e-3]
    else:
        coefficients = [-1.15282e-1, 6.63234e-2, -1.29163e-2, 8.77272e-4]

    calibration = CernoxCalibration(coefficients)
    temperature = calibration.temperature_from_resistance(resistance)

    st.write(f"Estimated temperature: {temperature:.2f} K")


if __name__ == "__main__":
    main()