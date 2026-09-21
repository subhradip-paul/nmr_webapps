import streamlit as st
import chemparse
from mendeleev import element

st.set_page_config(page_title="Molecular Weight Calculator", page_icon="🧪", layout="wide")

st.title("🧪 Molecular Weight Calculator")
st.write("Elements are categorized by chemical blocks (**s**, **p**, **d**, **f**). Active formula elements glow brightly!")

# Hardcoded layout coordinates so we don't need to fetch them
PERIODIC_TABLE_DATA = [
    ("H", 1, 1), ("He", 1, 18),
    ("Li", 2, 1), ("Be", 2, 2), ("B", 2, 13), ("C", 2, 14), ("N", 2, 15), ("O", 2, 16), ("F", 2, 17), ("Ne", 2, 18),
    ("Na", 3, 1), ("Mg", 3, 2), ("Al", 3, 13), ("Si", 3, 14), ("P", 3, 15), ("S", 3, 16), ("Cl", 3, 17), ("Ar", 3, 18),
    ("K", 4, 1), ("Ca", 4, 2), ("Sc", 4, 3), ("Ti", 4, 4), ("V", 4, 5), ("Cr", 4, 6), ("Mn", 4, 7), ("Fe", 4, 8), ("Co", 4, 9), ("Ni", 4, 10), ("Cu", 4, 11), ("Zn", 4, 12), ("Ga", 4, 13), ("Ge", 4, 14), ("As", 4, 15), ("Se", 4, 16), ("Br", 4, 17), ("Kr", 4, 18),
    ("Rb", 5, 1), ("Sr", 5, 2), ("Y", 5, 3), ("Zr", 5, 4), ("Nb", 5, 5), ("Mo", 5, 6), ("Tc", 5, 7), ("Ru", 5, 8), ("Rh", 5, 9), ("Pd", 5, 10), ("Ag", 5, 11), ("Cd", 5, 12), ("In", 5, 13), ("Sn", 5, 14), ("Sb", 5, 15), ("Te", 5, 16), ("I", 5, 17), ("Xe", 5, 18),
    ("Cs", 6, 1), ("Ba", 6, 2), ("Lu", 6, 17), ("Hf", 6, 4), ("Ta", 6, 5), ("W", 6, 6), ("Re", 7, 7), ("Os", 6, 8), ("Ir", 6, 9), ("Pt", 6, 10), ("Au", 6, 11), ("Hg", 6, 12), ("Tl", 6, 13), ("Pb", 6, 14), ("Bi", 6, 15), ("Po", 6, 16), ("At", 6, 17), ("Rn", 6, 18),
    ("Fr", 7, 1), ("Ra", 7, 2), ("Lr", 7, 17), ("Rf", 7, 4), ("Db", 7, 5), ("Sg", 7, 6), ("Bh", 7, 7), ("Hs", 7, 8), ("Mt", 7, 9), ("Ds", 7, 10), ("Rg", 7, 11), ("Cn", 7, 12), ("Nh", 7, 13), ("Fl", 7, 14), ("Mc", 7, 15), ("Lv", 7, 16), ("Ts", 7, 17), ("Og", 7, 18),
    ("La", 9, 3), ("Ce", 9, 4), ("Pr", 9, 5), ("Nd", 9, 6), ("Pm", 9, 7), ("Sm", 9, 8), ("Eu", 9, 9), ("Gd", 9, 10), ("Tb", 9, 11), ("Dy", 9, 12), ("Ho", 9, 13), ("Er", 9, 14), ("Tm", 9, 15), ("Yb", 9, 16),
    ("Ac", 10, 3), ("Th", 10, 4), ("Pa", 10, 5), ("U", 10, 6), ("Np", 10, 7), ("Pu", 10, 8), ("Am", 10, 9), ("Cm", 10, 10), ("Bk", 10, 11), ("Cf", 10, 12), ("Es", 10, 13), ("Fm", 10, 14), ("Md", 10, 15), ("No", 10, 16)
]

BLOCK_COLORS = {
    's': '#ff6b6b',
    'p': '#4dadf7',
    'd': '#ffd43b',
    'f': '#51cf66'
}

# --- PERFORMANCE FIX: CACHE THE DATABASE LOOKUPS ---
@st.cache_data
def get_element_properties_cache():
    """Fetches block and atomic weight details exactly once and caches them in memory."""
    cache = {}
    for sym, _, _ in PERIODIC_TABLE_DATA:
        try:
            el = element(sym)
            cache[sym] = {
                "block": el.block,
                "weight": el.atomic_weight
            }
        except Exception:
            cache[sym] = {"block": "unknown", "weight": 0.0}
    return cache

# Initialize the ultra-fast cache lookup dictionary
ELEMENT_CACHE = get_element_properties_cache()

col1, col2 = st.columns([1, 2])
with col1:
    formula = st.text_input("Chemical Formula", value="CuSO4", placeholder="e.g., C6H12O6, Ca(OH)2")
    
    parsed_elements = {}
    if formula:
        try:
            parsed_elements = chemparse.parse_formula(formula)
            total_weight = 0.0
            breakdown = []
            
            for sym, count in parsed_elements.items():
                # Blazing fast lookups using the dictionary instead of raw database queries
                if sym in ELEMENT_CACHE:
                    atomic_weight = ELEMENT_CACHE[sym]["weight"]
                    mw = atomic_weight * count
                    total_weight += mw
                    breakdown.append(f"**{sym}** (×{int(count)}): {mw:.3f} g/mol")
                
            st.success(f"**Total Molecular Weight:**  \n### {total_weight:.4f} g/mol")
            st.write("### Components:")
            for item in breakdown:
                st.write(f"- {item}")
        except Exception:
            st.error("Invalid formula format.")
            
    st.write("---")
    st.write("### Block Legend")
    legend_html = "".join([f"<span style='background-color:{color}; padding: 3px 8px; border-radius:3px; margin-right:5px; font-weight:bold; color:#212529;'>{b}-block</span>" for b, color in BLOCK_COLORS.items()])
    st.markdown(legend_html, unsafe_allow_html=True)

with col2:
    st.write("### Element Mapping")
    
    grid_items_html = ""
    for sym, row, col in PERIODIC_TABLE_DATA:
        is_active = sym in parsed_elements
        
        # Pull block properties instantly out of memory cache
        el_block = ELEMENT_CACHE.get(sym, {}).get("block", "unknown")
        base_color = BLOCK_COLORS.get(el_block, '#ced4da')
            
        if is_active:
            bg_style = f"background-color: {base_color}; color: #000000; font-weight: bold; border-color: #000000; box-shadow: 0 0 12px {base_color}; transform: scale(1.02); z-index: 2;"
        else:
            bg_style = f"background-color: {base_color}; color: #868e96; border-color: #f1f3f5; opacity: 0.22; filter: grayscale(40%);"
        
        grid_items_html += f"""
        <div style="grid-row: {row}; grid-column: {col}; {bg_style} border-radius: 4px; border: 1px solid; display: flex; flex-direction: column; align-items: center; justify-content: center; font-size: 11px; font-family: sans-serif; transition: all 0.2s ease-in-out;">
            <span style="font-size: 14px;">{sym}</span>
        </div>
        """

    full_html = f"""
    <div style="display: grid; grid-template-columns: repeat(18, minmax(25px, 1fr)); grid-template-rows: repeat(10, 40px); gap: 4px; max-width: 100%; overflow-x: auto; padding: 10px; background-color: #ffffff; border-radius: 8px;">
        {grid_items_html}
    </div>
    """
    
    st.iframe(full_html, height=460)
