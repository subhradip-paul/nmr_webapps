import streamlit as st
import numpy as np
import plotly.graph_objects as go


st.set_page_config()

log2 = np.log(2)

# ----- 1. Pseudo Voigt Function -----
def pvoigt_nmr(x, peak_height, peak_position, line_width, fraction):
    sigma = line_width / 2.0
    sigmag = sigma / np.sqrt(2 * log2)
    amplitude = (((peak_height) / 2) * np.pi * line_width) / (1 + fraction * (np.sqrt(np.pi * log2) - 1))
    g_nmr = (((1 - fraction) * amplitude) / (sigmag * np.sqrt(2 * np.pi))) * np.exp(-((x - peak_position) ** 2 / (2 * sigmag ** 2)))
    l_nmr = (1 / np.pi) * fraction * amplitude * sigma ** 1 / ((x - peak_position) ** 2 + sigma ** 2)
    return g_nmr + l_nmr

# ----- 2. DQ-SQ Generator -----
def gen2ddqsq(peak_posns, coup_matrix, cfreqf1, cfreqf2, swf1, swf2, fwhhf1, fwhhf2):
    xscale = np.linspace(cfreqf2 - swf2 / 2, cfreqf2 + swf2 / 2, 1000)
    yscale = np.linspace(-swf1 / 2 + cfreqf1, swf1 / 2 + cfreqf1, 1000)
    int_matrixv = np.zeros((xscale.size, yscale.size))
    folded_f1_posn = np.zeros(coup_matrix.shape)
    for idx1, ii in enumerate(coup_matrix):
        for idx2, jj in enumerate(ii):
            if jj > 0:
                peak_posnf2 = peak_posns[idx1]
                peak_posnf1 = np.mod(peak_posns[idx1] + peak_posns[idx2] + swf1 / 2 - cfreqf1, swf1) - swf1 / 2 + cfreqf1
                folded_f1_posn[idx1, idx2] = peak_posnf1
                voigtx = pvoigt_nmr(xscale, jj, peak_posnf2, fwhhf2, 0.2)
                voigty = pvoigt_nmr(yscale, jj, peak_posnf1, fwhhf1, 0.2)
                int_matrixv += np.outer(voigtx, voigty)
    folded_f1_posn = np.tril(folded_f1_posn)
    return int_matrixv, folded_f1_posn

if 'simulate_spectrum' not in st.session_state:
    st.session_state.simulate_spectrum = False

def onClickfn():
    st.session_state.simulate_spectrum = True


# ----- 3. User Inputs -----
st.title("2D DQ-SQ NMR Spectrum Generator")
st.divider()

st.write('''
The script helps you to generate a double-quantum, single-quantum correlation spectrum
and calculate where the peaks will be in case they are folded.
It should work with spectrum acquired with STATES-TPPI method.
''')

num_peaks = st.number_input("Number of Peaks", min_value=2, max_value=10, value=3, step=1)

st.subheader("Peak Names and Positions (ppm)")
peak_names, peak_posns = [], []
for i in range(num_peaks):
    cols = st.columns([1, 2])
    with cols[0]:
        name = st.text_input(f"Peak {i+1} Name", value=f"P{i+1}", key=f"name_{i}")
    with cols[1]:
        shift = st.number_input(f"Shift (ppm) for {name}", value=float(-110 + i*10), key=f"shift_{i}")
    peak_names.append(name)
    peak_posns.append(shift)

peak_posns = np.array(peak_posns)

st.subheader("Coupling Matrix (Z-Matrix)")

z_matrix = np.zeros((num_peaks, num_peaks))
for i in range(num_peaks):
    cols = st.columns(num_peaks)
    for j in range(i + 1):  # lower triangle including diagonal
        z_matrix[i][j] = cols[j].number_input(f"Z[{i+1},{j+1}]", min_value=0.0, max_value=1.0, value=0.5, step=0.1, key=f"z_{i}_{j}")
coup_matrix = z_matrix + z_matrix.T

st.subheader("Spectral Parameters")
cfreqf2 = st.number_input("Carrier frequency F2 (ppm)", value=np.mean(peak_posns))
cfreqf1 = st.number_input("Carrier frequency F1 (DQ, ppm)", value=2 * cfreqf2)
swf2 = st.number_input("Spectral width F2 (ppm)", value=np.max(peak_posns)-np.min(peak_posns)+50)
swf1 = st.number_input("Spectral width F1 (ppm)", value=2*swf2)

st.session_state.reverse = None

def on_toggle_change():
    if st.session_state.reverse:
        st.session_state.reverse = False
    else:
        st.session_state.reverse = True


# ----- 4. Generate and Plot -----
st.button("Simulate Spectrum", on_click = onClickfn)

if st.session_state.simulate_spectrum:

    line_width_f2 = st.number_input("Line width in F2 (ppm)", value=2.0)
    line_width_f1 = st.number_input("Line width in F1 (ppm)", value=2.0)

    int_matrixvfold, folded_f1_posn = gen2ddqsq(
        peak_posns, coup_matrix, cfreqf1, cfreqf2, swf1, swf2, line_width_f1, line_width_f2
    )

    xscale = np.linspace(cfreqf2 - swf2 / 2, cfreqf2 + swf2 / 2, 1000)
    yscale = np.linspace(-swf1 / 2 + cfreqf1, swf1 / 2 + cfreqf1, 1000)

    fig = go.Figure()
    color_contour = st.selectbox('Colorscale', (
    'aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance', 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl',
    'brbg', 'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl', 'darkmint', 'deep', 'delta', 'dense',
    'earth', 'edge', 'electric', 'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys', 'haline', 'hot', 'hsv',
    'ice', 'icefire', 'inferno', 'jet', 'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges', 'orrd',
    'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu',
    'pubugn', 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy', 'rdpu', 'rdylbu', 'rdylgn',
    'redor', 'reds', 'solar', 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose', 'tempo',
    'temps', 'thermal', 'tropic', 'turbid', 'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd'),
                              index=0, placeholder='Colour scheme for the contours'),

    reverse = st.toggle("Reverse colour map?", on_change=on_toggle_change())

    if reverse:
        color_scale = color_contour[0] + '_r'
        st.write(color_scale)
    else:
        color_scale = color_contour[0]

    # Add Contour plot
    fig.add_trace(go.Contour(
        z=int_matrixvfold.T,
        x=xscale,
        y=yscale,

        colorscale=color_scale,


        contours=dict(start=0.2, end=np.max(int_matrixvfold), size=(np.max(int_matrixvfold) - 0.2) / 50),
        # colorscale=[[0, 'white'], [0.5, 'mediumturquoise'], [1, 'darkturquoise']],
        line_smoothing=0.85,
        contours_coloring='heatmap',
        colorbar=dict(title="Intensity"),
        line=dict(width=0.01, color='limegreen'),
        showscale=True,
        hoverinfo='x+y+z'
    ))

    # Add diagonal lines
    for k in range(-2, 3):
        x0 = cfreqf2 - k * swf1 / 2
        x1 = cfreqf2 + swf2 / 2
        y0 = cfreqf1 + 2 * (x0 - cfreqf2)
        y1 = cfreqf1 + 2 * (x1 - cfreqf2)
        fig.add_trace(go.Scatter(
            x=[x0, x1],
            y=[y0, y1],
            mode="lines",
            line=dict(color='red', width=1, dash='dot'),
            showlegend=False,
        ))

    # Add horizontal connectivity lines
    for idx, shifts in np.ndenumerate(folded_f1_posn):
        if shifts != 0:
            xmin = peak_posns[idx[1]]
            xmax = peak_posns[idx[0]]
            yval = shifts
            fig.add_trace(go.Scatter(
                x=[xmin, xmax],
                y=[yval, yval],
                mode='lines',
                line=dict(color='black', width=0.25, dash='dot'),
                showlegend=False
            ))

    fig.update_layout(
        xaxis_title='SQ Dimension (ppm)',
        yaxis_title='DQ Dimension (ppm)',
        title='2D DQ-SQ NMR Spectrum (Plotly)',
        height=600,
        width=900,
        xaxis_range=[cfreqf2 + swf2 / 2, cfreqf2 - swf2 / 2],
        yaxis_range=[cfreqf1 + swf1 / 2, cfreqf1 - swf1 / 2],
        # xaxis=dict(autorange='reversed'),  # ppm axis flips
        # yaxis=dict(autorange='reversed'),
    )

    st.plotly_chart(fig, use_container_width=True)