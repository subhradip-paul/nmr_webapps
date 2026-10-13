import streamlit as st
import pandas as pd
import numpy as np
import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor

rdDepictor.SetPreferCoordGen(True)
from streamlit_ketcher import st_ketcher
from io import StringIO
import py3Dmol
from stmol import showmol, add_hover

script_dir = os.path.dirname(__file__)
csv_file = os.path.join(script_dir, '../dep/NMR_freq_table.csv')

nuctable = pd.read_csv(csv_file)
gyr_ratio_MHz_T = nuctable["GyrHz"]
name_nuc = nuctable["Name"]


def render3dmol(mol, style):
    mol_block = Chem.MolToMolBlock(mol)

    xyzview = py3Dmol.view()
    xyzview.addModel(mol_block, 'mol')

    if style == 'stick':
        xyzview.setStyle({'stick': {}})
    elif style == 'ball and stick':
        xyzview.setStyle({'sphere': {'scale': 0.3}, 'stick': {}})
    elif style == 'sphere':
        xyzview.setStyle({'sphere': {}})

    colour_bg = st.color_picker('Pick a background colour', '#9DEBF7')
    xyzview.setBackgroundColor(colour_bg)
    xyzview.zoomTo()
    add_hover(xyzview)
    showmol(xyzview, height=500, width=800)


def nuc_label(name, symbol):
    """Isotope label like '13C' for SLEEPY/SIMPSON. Assumes the 'Name' column
    looks like '1H' or '13C'; otherwise falls back to the element symbol."""
    m = re.fullmatch(r"\s*(\d+)\s*([A-Z][a-z]?)\s*", str(name))
    return f"{m.group(1)}{m.group(2)}" if m else str(symbol)


def get_gyration_inputs(atom_list):
    """UI for isotope selection. Returns (gyromagnetic ratios, isotope labels)."""
    gyr_atom = []
    labels = []

    st.subheader("Isotope Selection")

    for idx, symbol in enumerate(atom_list):
        matches = nuctable[nuctable["Symbol"] == symbol]

        if matches.shape[0] == 1:
            row = matches.iloc[0]
            gyr_atom.append(row["GyrHz"])
            labels.append(nuc_label(row["Name"], symbol))

        elif matches.shape[0] > 1:
            key = f"isotope_{idx}_{symbol}"
            names = matches["Name"].tolist()

            if st.session_state.get(key) not in names:
                st.session_state[key] = names[0]

            choice = st.selectbox(f"{symbol} (atom {idx+1})", names, key=key)

            gyr_atom.append(matches[matches["Name"] == choice]["GyrHz"].values[0])
            labels.append(nuc_label(choice, symbol))

        else:
            st.warning(f"{symbol} not found. Using 0.")
            gyr_atom.append(0.0)
            labels.append(symbol)

    return np.array(gyr_atom), labels


def dipolar_euler(ri, rj):
    """ZYZ Euler angles (deg) that take z onto the i->j internuclear vector.
    gamma is arbitrary for an axial tensor."""
    v = np.asarray(rj, float) - np.asarray(ri, float)
    v /= np.linalg.norm(v)
    beta = np.degrees(np.arccos(np.clip(v[2], -1, 1)))
    alpha = np.degrees(np.arctan2(v[1], v[0]))
    return np.array([alpha, beta, 0.0])


def xyz_to_dipolar_data(xyz_dataframe, gyr_atom):
    num_atoms = len(xyz_dataframe)

    nuc = xyz_dataframe['atom'].tolist()
    coord_xyz = xyz_dataframe[['x', 'y', 'z']].astype(float).to_numpy()

    pl = 6.62607015e-34

    results = []

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):

            r = np.linalg.norm(coord_xyz[i] - coord_xyz[j])

            gyr1 = gyr_atom[i] * 1e6
            gyr2 = gyr_atom[j] * 1e6

            dip = -1e-7 * (gyr1 * gyr2 * pl) / ((r * 1e-10) ** 3)

            euler = dipolar_euler(coord_xyz[i], coord_xyz[j])

            results.append([
                i + 1, nuc[i],
                j + 1, nuc[j],
                round(r, 2),
                round(dip, 2),
                *np.round(euler, 2)
            ])

    return pd.DataFrame(results, columns=[
        'i', 'Nuc i', 'j', 'Nuc j',
        'Distance (Å)', 'Dipolar Coupling (Hz)',
        'alpha', 'beta', 'gamma'
    ])


# ---------------------------------------------------------------------------
# Text output for SLEEPY and SIMPSON
# ---------------------------------------------------------------------------
def make_sleepy_text(df, nuclei):
    """SLEEPY: 0-based indices, delta = 2*b (b in Hz)."""
    lines = [f"nuclei = {list(nuclei)}", ""]
    for _, r in df.iterrows():
        lines.append(
            f"ex.set_inter('dipole', i0={int(r['i']) - 1}, i1={int(r['j']) - 1}, "
            f"delta={2 * r['Dipolar Coupling (Hz)']:.2f}, "
            f"euler_d=[{r['alpha']:.2f}, {r['beta']:.2f}, {r['gamma']:.2f}])"
        )
    return "\n".join(lines)


def make_simpson_text(df, nuclei):
    """SIMPSON: 1-based indices, dipole i j b(Hz) alpha beta gamma."""
    channels = list(dict.fromkeys(nuclei))
    lines = ["spinsys {",
             f"  channels {' '.join(channels)}",
             f"  nuclei {' '.join(nuclei)}"]
    for _, r in df.iterrows():
        lines.append(
            f"  dipole {int(r['i'])} {int(r['j'])} "
            f"{r['Dipolar Coupling (Hz)']:.2f} "
            f"{r['alpha']:.2f} {r['beta']:.2f} {r['gamma']:.2f}"
        )
    lines.append("}")
    return "\n".join(lines)


def show_dipolar_output(df, nuclei):
    """Dipolar table with row selection, plus copy-ready SLEEPY/SIMPSON text.
    With no rows selected, all pairs are exported."""
    st.write("### Dipolar couplings")
    st.caption("Select rows to export only those pairs (none selected = all pairs).")
    event = st.dataframe(
        df,
        hide_index=True,
        width='content',
        on_select="rerun",
        selection_mode="multi-row",
    )
    rows = event.selection.rows
    export_df = df.iloc[rows] if rows else df

    st.subheader("SLEEPY input")
    st.code(make_sleepy_text(export_df, nuclei), language="python")

    st.subheader("SIMPSON input")
    st.code(make_simpson_text(export_df, nuclei), language="text")


def dist2dipole(choice_nuc1, choice_nuc2, distance):
    pl = 6.62607015e-34
    nuc1idx = name_nuc[name_nuc == choice_nuc1].index
    nuc2idx = name_nuc[name_nuc == choice_nuc2].index

    gyr1 = gyr_ratio_MHz_T[nuc1idx[0]] * 1e6
    gyr2 = gyr_ratio_MHz_T[nuc2idx[0]] * 1e6

    dip = -1e-7 * (gyr1 * gyr2 * pl) / ((distance * 1e-10) ** 3)
    return dip


def dipole2dist(choice_nuc1, choice_nuc2, dipole):
    pl = 6.62607015e-34
    nuc1idx = name_nuc[name_nuc == choice_nuc1].index
    nuc2idx = name_nuc[name_nuc == choice_nuc2].index

    gyr1 = gyr_ratio_MHz_T[nuc1idx[0]] * 1e6
    gyr2 = gyr_ratio_MHz_T[nuc2idx[0]] * 1e6

    dist = 1e10 * ((1e-7 * abs((gyr1 * gyr2 * pl)) / dipole) ** (1 / 3))
    return dist


st.title(r"Dipolar Coupling Calculator")
st.divider()

choice_of_calculation = st.radio(
    ":rainbow[**What do you want to calculate?**]",
    [r"Dipole :arrows_counterclockwise: Distance", r'Dipolar Couplings from Structure'],
    index=None)

if choice_of_calculation == "Dipole :arrows_counterclockwise: Distance":
    choice_nuc1 = st.selectbox('Nucleus 1:', nuctable.Name.values.tolist(), index=2)
    choice_nuc2 = st.selectbox('Nucleus 2:', nuctable.Name.values.tolist(), index=11)
    choice_calc = st.radio('Method', ('Distance to Dipole', 'Dipole to Distance'))

    if choice_calc == 'Distance to Dipole':
        distance = st.number_input('Enter distance in Angstrom: ', min_value=0.5, value=1.07)
        coup = dist2dipole(choice_nuc1, choice_nuc2, distance)
        st.divider()
        st.write("The dipolar coupling = " + str(np.abs(np.round(coup, 2))) + ' Hz')
        st.write("The dipolar coupling = " + str(np.abs(np.round(coup / 1e3, 2))) + ' kHz')
    else:
        dipole = st.number_input('Enter dipolar coupling in Hz: ', min_value=0.5)
        d = dipole2dist(choice_nuc1, choice_nuc2, dipole)
        st.divider()
        st.latex(r"\text{Distance  = }" + str(np.abs(np.round(d, 2))) + r'\ {\AA}')

elif choice_of_calculation == 'Dipolar Couplings from Structure':

    options = st.selectbox(
        "Draw a structure or manually enter *xyz* coordinates",
        options=["Draw", "Manual"],
        index=0
    )

    if options == "Manual":

        xyz_lines = st.file_uploader("Upload an xyz coordinates file", type="xyz")

        if xyz_lines is not None:

            df = pd.read_csv(
                xyz_lines,
                skiprows=2,
                sep=r"\s+",
                names=["atom", "x", "y", "z"],
                dtype={"atom": str, "x": float, "y": float, "z": float}
            )

            st.write("### Input structure")
            st.dataframe(df)

            gyr_atom, nuc_labels = get_gyration_inputs(df["atom"].tolist())

            df_xyz_dipole = xyz_to_dipolar_data(xyz_dataframe=df, gyr_atom=gyr_atom)

            show_dipolar_output(df_xyz_dipole, nuc_labels)

    else:

        smiles = st_ketcher("cc")

        if smiles is not None:

            mol = Chem.MolFromSmiles(smiles)

            if mol is not None:

                remove_1H = st.checkbox("Remove 1H")

                mol3d = Chem.AddHs(mol)

                params = AllChem.ETKDGv3()
                params.randomSeed = 0xf00d

                status = AllChem.EmbedMolecule(mol3d, params)

                if status != 0:
                    st.error("Could not generate a 3D structure.")
                    st.stop()

                AllChem.MMFFOptimizeMolecule(mol3d)

                if remove_1H:
                    mol3d = Chem.RemoveAllHs(mol3d)

                style = st.selectbox("Select style:", ["stick", "ball and stick", "sphere"])

                render3dmol(mol3d, style)

                xyz_string = Chem.MolToXYZBlock(mol3d)
                xyz_lines = xyz_string.strip().splitlines()

                df = pd.read_csv(
                    StringIO("\n".join(xyz_lines[2:])),
                    sep=r"\s+",
                    names=["atom", "x", "y", "z"]
                )

                # A checkbox instead of st.button: a button resets on every
                # rerun, which would hide the output whenever a selectbox or
                # table row below it is changed.
                if st.checkbox("Calculate"):

                    st.write("### Input structure")
                    st.dataframe(df)

                    gyr_atom, nuc_labels = get_gyration_inputs(df["atom"].tolist())

                    df_xyz_dipole = xyz_to_dipolar_data(xyz_dataframe=df, gyr_atom=gyr_atom)

                    show_dipolar_output(df_xyz_dipole, nuc_labels)