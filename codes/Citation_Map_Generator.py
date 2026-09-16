"""
Streamlit app: Zotero collection -> OpenAlex forward citation tree

Pipeline:
    1. Pull publication items from a Zotero collection via the Zotero Web API.
    2. Resolve each item's DOI to an OpenAlex work (cached).
    3. Let the user pick a publication from a selectbox.
    4. Build a forward citation tree (papers citing it, then papers citing
       those) and render it as an interactive Plotly graph. Clicking a node
       shows a link to open that paper.

Setup:
    pip install streamlit pyzotero requests plotly networkx --break-system-packages

    Get your Zotero credentials:
      - Library ID:  zotero.org -> Settings -> Feeds/API -> "Your userID for use in API calls"
      - API key:     zotero.org/settings/keys -> Create new private key (read-only is enough)
      - Collection key: open the collection in the Zotero web library; it's in the URL
                        (or leave ZOTERO_COLLECTION_KEY unset to pull the whole library)

    Store credentials in codes/.streamlit/secrets.toml (never hardcode them, never commit):
        ZOTERO_LIBRARY_ID = "1234567"
        ZOTERO_LIBRARY_TYPE = "user"   # or "group" for a shared lab library
        ZOTERO_API_KEY = "xxxxxxxxxxxxxxxxxxxx"
        ZOTERO_COLLECTION_KEY = "ABCD1234"   # optional
        OPENALEX_EMAIL = "your@email.com"

Run:
    streamlit run Citation_Map_Generator.py
"""

import streamlit as st
from pyzotero import zotero
import requests
import plotly.graph_objects as go
import networkx as nx

OPENALEX_BASE = "https://api.openalex.org"

# Depth-based styling: 0 = seed, 1 = direct citations, 2 = citations of citations, ...
DEPTH_COLORS = {0: "#d62728", 1: "#2ca02c", 2: "#1f77b4", 3: "#9467bd"}
DEPTH_LABELS = {0: "Selected paper", 1: "Cites it directly", 2: "Cites a citing paper", 3: "3rd generation"}


# ---------- Zotero ----------

@st.cache_data(ttl=3600)  # re-check Zotero once an hour
def fetch_zotero_items():
    """Pull items (with DOIs) from the configured Zotero library/collection."""
    zot = zotero.Zotero(
        st.secrets["ZOTERO_LIBRARY_ID"],
        st.secrets["ZOTERO_LIBRARY_TYPE"],
        st.secrets["ZOTERO_API_KEY"],
    )

    collection_key = st.secrets.get("ZOTERO_COLLECTION_KEY")
    if collection_key:
        items = zot.everything(zot.collection_items(collection_key, itemType="-attachment"))
    else:
        items = zot.everything(zot.top(itemType="-attachment"))

    pubs = []
    for item in items:
        data = item.get("data", {})
        doi = data.get("DOI")
        if not doi:
            continue  # skip items without a DOI — can't resolve to OpenAlex reliably
        pubs.append({
            "title": data.get("title", "Untitled"),
            "doi": doi,
            "year": data.get("date", "")[:4],
        })
    return pubs


# ---------- OpenAlex ----------

@st.cache_data(ttl=86400)  # citation data barely changes day to day
def resolve_openalex_work(doi):
    """Look up a single DOI on OpenAlex. Returns None if not found."""
    email = st.secrets.get("OPENALEX_EMAIL", "")
    try:
        r = requests.get(
            f"{OPENALEX_BASE}/works/doi:{doi}",
            params={"mailto": email},
            timeout=20,
        )
        if r.status_code != 200:
            return None
        return r.json()
    except requests.RequestException:
        return None


def _short_label(work):
    authors = work.get("authorships", [])
    first_author = authors[0]["author"]["display_name"] if authors else "Unknown"
    last_name = first_author.split()[-1] if first_author else "?"
    year = work.get("publication_year", "?")
    return f"{last_name} {year}"


def _fetch_citing(work_id, limit, email):
    r = requests.get(
        f"{OPENALEX_BASE}/works",
        params={"filter": f"cites:{work_id}", "per_page": limit, "mailto": email},
        timeout=20,
    )
    if r.status_code != 200:
        return []
    return r.json().get("results", [])


@st.cache_data(ttl=86400)
def build_forward_citation_graph(seed_id, max_depth=2, max_per_level_by_depth=None, max_total_nodes=60):
    """
    Breadth-first expansion of the forward citation tree: seed -> papers
    citing it -> papers citing those, out to max_depth generations.

    max_per_level_by_depth: dict mapping depth -> how many citing papers to
        fetch per node at that depth, e.g. {1: 10, 2: 2}. Depths not listed
        fall back to 5.
    max_total_nodes: hard cap on total nodes across the whole tree.

    Returns (nodes_dict, edges_list). Each node has label/title/doi/depth.
    Edges point citing_paper -> cited_paper, matching the citation direction.
    """
    max_per_level_by_depth = max_per_level_by_depth or {}
    email = st.secrets.get("OPENALEX_EMAIL", "")
    seed = requests.get(f"{OPENALEX_BASE}/works/{seed_id}", params={"mailto": email}).json()

    nodes = {seed_id: {"label": _short_label(seed), "title": seed["display_name"],
                        "doi": seed.get("doi"), "depth": 0}}
    edges = []
    frontier = [seed_id]

    for depth in range(1, max_depth + 1):
        limit = max_per_level_by_depth.get(depth, 5)
        next_frontier = []
        for parent_id in frontier:
            if len(nodes) >= max_total_nodes:
                break
            citing_works = _fetch_citing(parent_id, limit, email)
            for w in citing_works:
                wid = w["id"].split("/")[-1]
                edges.append((wid, parent_id))  # citing paper -> paper it cites
                if wid not in nodes:
                    nodes[wid] = {"label": _short_label(w), "title": w["display_name"],
                                  "doi": w.get("doi"), "depth": depth}
                    next_frontier.append(wid)
                if len(nodes) >= max_total_nodes:
                    break
            if len(nodes) >= max_total_nodes:
                break
        frontier = next_frontier
        if not frontier:
            break  # no more citing papers found at this depth

    return nodes, edges


# ---------- Plotly rendering ----------

def _node_url(node_id, nodes_dict):
    doi = nodes_dict[node_id].get("doi")
    return doi if doi else f"https://openalex.org/{node_id}"


def render_citation_graph(nodes_dict, edges_list):
    """nodes_dict: {id: {"label":..., "title":..., "doi":..., "depth":...}}
       edges_list: [(source_id, target_id), ...]

    Renders the graph and, if the user clicks a node, shows a link to open
    that paper (DOI if available, else its OpenAlex page).
    """
    G = nx.DiGraph()
    G.add_nodes_from(nodes_dict.keys())
    G.add_edges_from(edges_list)
    pos = nx.spring_layout(G, seed=42, k=0.7)

    edge_x, edge_y = [], []
    for u, v in edges_list:
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=1, color="#999999"),
                             hoverinfo="none", mode="lines")

    node_traces = []
    depths_present = sorted(set(n["depth"] for n in nodes_dict.values()))
    for depth in depths_present:
        ids = [nid for nid, n in nodes_dict.items() if n["depth"] == depth]
        node_traces.append(go.Scatter(
            x=[pos[nid][0] for nid in ids],
            y=[pos[nid][1] for nid in ids],
            mode="markers+text",
            text=[nodes_dict[nid]["label"] for nid in ids],
            textposition="top center",
            hovertext=[nodes_dict[nid]["title"] + "<br>Click to open" for nid in ids],
            hoverinfo="text",
            customdata=ids,
            marker=dict(size=22 if depth == 0 else 12,
                        color=DEPTH_COLORS.get(depth, "#7f7f7f"),
                        line=dict(width=1, color="white")),
            name=DEPTH_LABELS.get(depth, f"Generation {depth}"),
        ))

    fig = go.Figure(data=[edge_trace, *node_traces])
    fig.update_layout(
        showlegend=True,
        hovermode="closest",
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=650,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
    )

    event = st.plotly_chart(fig, width="stretch", on_select="rerun",
                             selection_mode="points", key="citation_graph")

    points = event.get("selection", {}).get("points", []) if event else []
    if points:
        clicked_id = points[0]["customdata"]
        clicked = nodes_dict[clicked_id]
        st.link_button(f"Open: {clicked['title']}", _node_url(clicked_id, nodes_dict))
    else:
        st.caption("Click a node to get a link to that paper.")


# ---------- Streamlit UI ----------

st.set_page_config(page_title="Lab Citation Map", layout="wide")
st.title("Lab Publications — Citation Map")

with st.spinner("Loading publications from Zotero..."):
    pubs = fetch_zotero_items()

if not pubs:
    st.error("No DOI-bearing items found. Check your Zotero credentials/collection key in secrets.toml.")
    st.stop()

pub_labels = [f"{p['title']} ({p['year']})" for p in pubs]
selected_idx = st.selectbox("Select a publication", range(len(pubs)), format_func=lambda i: pub_labels[i])
selected_pub = pubs[selected_idx]

with st.spinner("Resolving on OpenAlex..."):
    work = resolve_openalex_work(selected_pub["doi"])

if not work:
    st.warning("Couldn't find this publication on OpenAlex (DOI may be too recent or unindexed).")
    st.stop()

seed_id = work["id"].split("/")[-1]
st.header(f"Citations on OpenAlex: {work.get('cited_by_count', 'n/a')}", width="stretch") 

direct_citations_n = st.number_input(
    "How many direct citations to show?",
    min_value=1, max_value=50, value=5, step=1,
    help="Papers citing this one directly. Citations-of-citations are capped at 2 per paper to keep the graph readable.",
)

with st.spinner("Building citation tree..."):
    nodes_dict, edges_list = build_forward_citation_graph(
        seed_id,
        max_depth=2,
        max_per_level_by_depth={1: direct_citations_n, 2: 2},
    )

# st.markdown("🔴 Selected paper &nbsp;&nbsp; 🟢 Cites it directly &nbsp;&nbsp; 🔵 Cites a citing paper")
render_citation_graph(nodes_dict, edges_list)