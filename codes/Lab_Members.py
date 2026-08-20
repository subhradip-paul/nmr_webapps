import os
import pandas as pd
import plotly.express as px
import pycountry
import streamlit as st

st.set_page_config(page_title="Lab Origins Map", layout="wide")

# ---------------------------------------------------------------------------
# 1. DATA SOURCE — Google Sheet published as CSV
# ---------------------------------------------------------------------------
# In Google Sheets: File > Share > Publish to web > select the sheet tab > CSV
# You'll get a URL like:
#   https://docs.google.com/spreadsheets/d/e/2PACX-XXXXXXX/pub?output=csv
# Set it as a Heroku config var (recommended) so you never hardcode it:
#   heroku config:set LAB_ROSTER_SHEET_URL="https://docs.google.com/.../pub?output=csv"

SHEET_CSV_URL = os.environ.get(
    "LAB_ROSTER_SHEET_URL",
    "https://docs.google.com/spreadsheets/d/e/2PACX-1vT-Gu_f-HkoUyUAXkSQ_DUZzplc39QdapFSHKtqvZdfbnnWEpkLLBUQ2qr2y1yLiv4V5pe24HwB1-mW/pub?gid=325625685&single=true&output=csv",
)

@st.cache_data(ttl=300)  # re-fetch at most every 5 minutes
def load_roster(url: str) -> pd.DataFrame:
    df = pd.read_csv(url)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


try:
    people = load_roster(SHEET_CSV_URL)
except Exception as e:
    st.error(f"Could not load roster from Google Sheet: {e}")
    st.stop()

if "country" not in people.columns:
    st.error("Sheet must have a 'country' column.")
    st.stop()

# # ---------------------------------------------------------------------------
# # 2. Optional filter by status (current vs alumni), if that column exists
# # ---------------------------------------------------------------------------

# if "status" in people.columns:
#     statuses = sorted(people["status"].dropna().unique().tolist())
#     chosen = st.multiselect("Show", statuses, default=statuses)
#     people = people[people["status"].isin(chosen)]
 
# ---------------------------------------------------------------------------
# 3. Aggregate counts per country (total + breakdown by status, if present)
# ---------------------------------------------------------------------------
counts = people.groupby("country").size().reset_index(name="num_people")
 
if "status" in people.columns:
    # Pivot to one column per status value, e.g. "current", "alumni"
    status_counts = (
        people.groupby(["country", "status"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    counts = counts.merge(status_counts, on="country", how="left")
    status_cols = [c for c in status_counts.columns if c != "country"]
else:
    status_cols = []

# # ---------------------------------------------------------------------------
# # 3. Aggregate counts per country
# # ---------------------------------------------------------------------------
# counts = people.groupby("country").size().reset_index(name="num_people")


# ---------------------------------------------------------------------------
# 4. Convert to ISO-3 codes for robust matching
# ---------------------------------------------------------------------------
def to_iso3(name):
    try:
        return pycountry.countries.lookup(name).alpha_3
    except LookupError:
        return None


counts["iso3"] = counts["country"].apply(to_iso3)
missing = counts[counts["iso3"].isna()]
if not missing.empty:
    st.warning(f"Couldn't match these country names: {missing['country'].tolist()}")

# ---------------------------------------------------------------------------
# 5. Build the choropleth
# ---------------------------------------------------------------------------
fig = px.choropleth(
    counts,
    locations="iso3",
    locationmode="ISO-3",
    color="num_people",
    hover_name="country",
    hover_data={
        "iso3": False,  # hide raw code from tooltip
        "num_people": True,
        **{col: True for col in status_cols},
    },
    color_continuous_scale="sunsetdark",
    projection="natural earth",
    labels={"num_people": "Total", **{col: col.title() for col in status_cols}},
)
 
fig.update_layout(
    margin=dict(l=0, r=0, t=30, b=0),
    coloraxis_colorbar=dict(title="People"),
    geo=dict(
        showframe=False,
        showcoastlines=True,
        bgcolor="rgba(0,0,0,0)",
        # Crop out Antarctica (and empty polar whitespace) by limiting
        # the visible latitude range rather than trying to hide one
        # country's shading specifically.
        lataxis=dict(range=[-58, 85]),
    ),
)
 
# ---------------------------------------------------------------------------
# 6. Layout
# ---------------------------------------------------------------------------
st.title("Where our lab members are from")
st.plotly_chart(fig, use_container_width=True)

with st.expander("Show raw data"):
    st.dataframe(people.sort_values("country"))

st.caption("Some people have rejoined the lab later with different statuses.")