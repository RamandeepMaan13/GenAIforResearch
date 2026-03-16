"""
app.py — AI Biological Discovery Assistant
Streamlit frontend for gene neglect scoring and hypothesis generation.
Ollama / Mistral backend only.
"""
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from collections import defaultdict

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Gene Discovery Assistant",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }
.stApp { background: #0d1117; color: #e6edf3; }
h1, h2, h3 { font-family: 'IBM Plex Mono', monospace !important; color: #58a6ff !important; }

.metric-card {
    background: #161b22; border: 1px solid #30363d;
    border-radius: 8px; padding: 16px 20px; margin: 4px 0;
}
.metric-label {
    font-size: 11px; color: #8b949e; text-transform: uppercase;
    letter-spacing: 0.08em; font-family: 'IBM Plex Mono', monospace;
}
.metric-value {
    font-size: 28px; font-weight: 600; color: #58a6ff;
    font-family: 'IBM Plex Mono', monospace; line-height: 1.2;
}

/* Big neglect score card */
.neglect-card {
    background: linear-gradient(135deg, #1a0e00, #2d1800);
    border: 1px solid #f97316; border-radius: 10px;
    padding: 20px 24px; margin: 6px 0; text-align: center;
}
.neglect-card .nlabel {
    font-size: 10px; color: #f97316; text-transform: uppercase;
    letter-spacing: 0.14em; font-family: 'IBM Plex Mono', monospace;
}
.neglect-card .nvalue {
    font-size: 40px; font-weight: 600; color: #fb923c;
    font-family: 'IBM Plex Mono', monospace; line-height: 1.1;
}
.neglect-card .nsub {
    font-size: 11px; color: #92400e;
    font-family: 'IBM Plex Mono', monospace; margin-top: 4px;
}

/* Signal breakdown rows */
.sig-row {
    display: flex; align-items: center; gap: 10px;
    margin: 6px 0; font-family: 'IBM Plex Mono', monospace; font-size: 12px;
}
.sig-label { color: #8b949e; width: 100px; flex-shrink: 0; }
.sig-bar-bg { flex: 1; background: #21262d; border-radius: 3px; height: 5px; }
.sig-bar-fill { height: 5px; border-radius: 3px; }
.sig-val { color: #c9d1d9; width: 50px; text-align: right; flex-shrink: 0; }

.scan-log {
    font-family: 'IBM Plex Mono', monospace; font-size: 11px;
    color: #8b949e; background: #0d1117; border: 1px solid #21262d;
    border-radius: 6px; padding: 12px; height: 200px;
    overflow-y: auto; line-height: 1.8;
}
.section-header {
    font-family: 'IBM Plex Mono', monospace; font-size: 10px;
    text-transform: uppercase; letter-spacing: 0.12em; color: #8b949e;
    border-bottom: 1px solid #21262d; padding-bottom: 6px; margin: 16px 0 8px 0;
}
.hypothesis-section {
    background: #161b22; border-left: 3px solid #58a6ff;
    padding: 10px 14px; margin: 8px 0; border-radius: 0 6px 6px 0;
    font-size: 14px; line-height: 1.6; color: #c9d1d9;
}
.hypothesis-section.gap        { border-color: #f97316; }
.hypothesis-section.experiment { border-color: #22c55e; }
.hypothesis-section.disease    { border-color: #a855f7; }
.hypothesis-section.drug       { border-color: #eab308; }

.rating-bar-bg {
    background: #21262d; border-radius: 4px; height: 8px;
    width: 100%; margin: 4px 0 10px 0;
}
.rating-bar-fill {
    height: 8px; border-radius: 4px;
    background: linear-gradient(90deg, #58a6ff, #3fb950);
}
.pathway-tag {
    display: inline-block; background: #1c2a1c; border: 1px solid #2d5a2d;
    color: #4ade80; border-radius: 3px; padding: 1px 6px;
    font-size: 11px; font-family: 'IBM Plex Mono', monospace; margin: 2px 1px;
}
.theme-tag {
    display: inline-block; background: #1a1a2e; border: 1px solid #3b3b8a;
    color: #818cf8; border-radius: 3px; padding: 1px 6px;
    font-size: 10px; font-family: 'IBM Plex Mono', monospace; margin: 1px;
}
</style>
""", unsafe_allow_html=True)

# ── Imports ───────────────────────────────────────────────────────────────────
from global_pathway_scan import scan_pathways
from gene_score_engine import compute_gene_score
from ai_engine import generate_hypothesis, parse_hypothesis

# ── Biological themes ─────────────────────────────────────────────────────────
# Keywords match against KEGG pathway names (e.g. "PI3K-Akt signaling pathway",
# "Fatty acid degradation", "Alzheimer disease") — keep them short and lowercase.
BIOLOGICAL_THEMES = {
    "All pathways":             None,
    "Cancer":                   ["cancer", "carcinoma", "melanoma", "glioma", "leukemia",
                                  "hepatocellular", "colorectal", "breast", "prostate", "lung"],
    "Metabolism":               ["metabol", "glycolysis", "fatty acid", "oxidative phosphorylation",
                                  "TCA cycle", "biosynthesis", "lipid", "amino acid", "carbon",
                                  "gluconeogenesis", "pentose", "butanoate", "propanoate"],
    "Aging & Longevity":        ["longevity", "FoxO", "mTOR", "AMPK", "PI3K-Akt",
                                  "cellular senescence", "autophagy", "sirtu"],
    "Immune System":            ["immune", "T cell", "B cell", "NK cell", "cytokine",
                                  "interleukin", "toll-like", "NF-kappa", "complement",
                                  "hematopoietic", "leukocyte"],
    "Neurodegeneration":        ["Alzheimer", "Parkinson", "Huntington", "prion",
                                  "neurodegeneration", "ALS", "spinocerebellar"],
    "Cardiovascular":           ["cardiac", "heart", "vascular", "hypertension",
                                  "atherosclerosis", "arrhythmia", "dilated cardiomyopathy",
                                  "hypertrophic cardiomyopathy"],
    "Cell Cycle & Apoptosis":   ["cell cycle", "apoptosis", "p53", "DNA repair",
                                  "cellular senescence", "ferroptosis", "necroptosis"],
    "Diabetes & Obesity":       ["insulin", "diabetes", "adipogenesis", "PPAR",
                                  "glucagon", "AGE-RAGE"],
    "Inflammation":             ["NF-kappa", "JAK-STAT", "TNF", "IL-17",
                                  "NOD-like", "RIG-I", "cytosolic DNA"],
    "Stem Cells & Development": ["stem cell", "Wnt", "Notch", "Hedgehog",
                                  "hippo", "TGF-beta", "ErbB"],
    "Epigenetics":              ["chromatin", "histone"],
    "Infectious Disease":       ["infection", "tuberculosis", "malaria", "influenza",
                                  "HIV", "hepatitis", "coronavirus", "Staphylococcus",
                                  "Salmonella", "pathogenic"],
    "Custom keyword":           "custom",
}

# ── Session state ─────────────────────────────────────────────────────────────
for key, default in {
    "scan_results": None,
    "scan_log":     [],
    "scan_done":    False,
    "ai_outputs":   {},
}.items():
    if key not in st.session_state:
        st.session_state[key] = default

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🧬 Gene Discovery")
    st.markdown("---")

    mode = st.radio(
        "Mode",
        ["Discovery Scan", "Single Gene Analysis"],
        label_visibility="collapsed"
    )

    st.markdown("### Biological Theme")
    theme_choice = st.selectbox(
        "Filter pathways by theme",
        options=list(BIOLOGICAL_THEMES.keys()),
        index=0,
        help="Narrows the KEGG pathway scan to a biological area of interest."
    )
    custom_keyword = ""
    if theme_choice == "Custom keyword":
        custom_keyword = st.text_input("Enter keyword", placeholder="e.g. autophagy, ferroptosis")
    elif theme_choice != "All pathways":
        kws = BIOLOGICAL_THEMES[theme_choice]
        st.markdown(
            "<div style='margin-top:4px'>"
            + " ".join(f'<span class="theme-tag">{k}</span>' for k in kws[:6])
            + "</div>",
            unsafe_allow_html=True
        )

    st.markdown("### Scan Settings")
    max_pathways = st.slider("Max pathways to scan", 5, 330, 20,
        help="330 = all human KEGG pathways. Start with 10–20 for speed.")
    top_n = st.slider("Top genes to show", 5, 30, 10)
    max_workers = st.slider("Parallel workers", 1, 10, 5,
        help="More workers = faster scoring but higher API load.")

    st.markdown("### Gene Source")
    scan_type = st.selectbox("Input type", ["KEGG pathway scan", "Custom gene list"])
    custom_genes_input = ""
    if scan_type == "Custom gene list":
        custom_genes_input = st.text_area(
            "Genes (comma separated)",
            placeholder="TP53, BRCA1, AKT1, FOXO3",
            height=80
        )

    st.markdown(
        "<div style='font-size:11px;color:#8b949e;margin-top:12px'>"
        "🤖 AI: Ollama / Mistral (local)</div>",
        unsafe_allow_html=True
    )
    st.markdown("---")
    if st.button("🗑️ Clear results", use_container_width=True):
        for k in ["scan_results", "scan_log", "scan_done", "ai_outputs"]:
            st.session_state[k] = [] if k in ("scan_log",) else (
                {} if k == "ai_outputs" else (False if k == "scan_done" else None)
            )
        st.rerun()

# ── Theme filter ──────────────────────────────────────────────────────────────
def pathway_matches_theme(pathway_name: str) -> bool:
    if theme_choice == "All pathways":
        return True
    if theme_choice == "Custom keyword":
        kw = custom_keyword.strip().lower()
        return bool(kw) and kw in pathway_name.lower()
    return any(k.lower() in pathway_name.lower() for k in BIOLOGICAL_THEMES[theme_choice])

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("# Gene Discovery Assistant")
st.markdown(
    "Score genes by **neglect** — high biological importance relative to research attention "
    "— then generate AI-powered hypotheses with Mistral."
)

# ═════════════════════════════════════════════════════════════════════════════
# DISCOVERY SCAN MODE
# ═════════════════════════════════════════════════════════════════════════════
if mode == "Discovery Scan":

    run_scan = st.button("▶  Run Discovery Scan", type="primary")

    if run_scan:
        st.session_state.scan_results = None
        st.session_state.scan_log     = []
        st.session_state.scan_done    = False
        st.session_state.ai_outputs   = {}

        log_placeholder = st.empty()
        progress_bar    = st.progress(0)
        status_text     = st.empty()
        log_lines       = []

        def progress_callback(current, total, message):
            log_lines.append(f"[{current:>4}/{total}] {message}")
            if len(log_lines) > 100:
                log_lines.pop(0)
            progress_bar.progress(min(int(current / max(total, 1) * 100), 100))
            status_text.markdown(f"`{message}`")
            log_html = "<br>".join(
                f"<span style='color:#3fb950'>›</span> {l}" for l in log_lines[-14:]
            )
            log_placeholder.markdown(
                f'<div class="scan-log">{log_html}</div>', unsafe_allow_html=True
            )

        # Custom gene list
        if scan_type == "Custom gene list":
            genes = [g.strip().upper() for g in custom_genes_input.split(",") if g.strip()]
            if not genes:
                st.error("Enter at least one gene name.")
                st.stop()
            results = []
            for i, gene in enumerate(genes):
                progress_callback(i + 1, len(genes), f"Scoring {gene}...")
                data = compute_gene_score(gene)
                if data:
                    data.update({"pathway_ids": ["Custom"], "pathway_names": ["Custom"], "pathway_count": 1})
                    results.append(data)
            st.session_state.scan_results = sorted(results, key=lambda x: x["score"], reverse=True)

        # KEGG scan
        else:
            # Pass theme filter into scan so KEGG pathways are filtered BEFORE
            # scoring — avoids scanning hundreds of irrelevant pathways
            theme_filter = pathway_matches_theme if theme_choice != "All pathways" else None
            all_results = scan_pathways(
                max_pathways=max_pathways,
                max_workers=max_workers,
                progress_callback=progress_callback,
                pathway_filter=theme_filter
            )
            st.session_state.scan_results = all_results

        progress_bar.progress(100)
        n = len(st.session_state.scan_results or [])
        if n == 0 and theme_choice != "All pathways":
            status_text.markdown(
                f"⚠️ No genes found for theme **{theme_choice}**. "
                f"Try a broader theme or use 'All pathways'."
            )
        else:
            status_text.markdown(
                f"✅ Scan complete — **{n} genes scored**"
                + (f" · theme: *{theme_choice}*" if theme_choice != "All pathways" else "")
            )
        st.session_state.scan_done = True

    # ── Results ───────────────────────────────────────────────────────────────
    if st.session_state.scan_results:
        results = st.session_state.scan_results
        display = results[:top_n]
        df      = pd.DataFrame(display)

        # Percentile rank across full result set
        all_scores = [r["score"] for r in results]
        def percentile_rank(score):
            below = sum(1 for s in all_scores if s < score)
            return round(below / max(len(all_scores), 1) * 100)

        # Signal bar maxima — normalize relative to the displayed gene set
        sig_max_c = max((r["centrality"] for r in display), default=1) or 1
        sig_max_e = max((r["entropy"]    for r in display), default=1) or 1
        sig_max_p = max((r["papers"]     for r in display), default=1) or 1
        sig_max_t = max((r["trend"]      for r in display), default=1) or 1

        st.markdown("---")

        # Summary metrics
        c1, c2, c3, c4, c5 = st.columns(5)
        with c1:
            st.markdown(f'<div class="metric-card"><div class="metric-label">Genes scored</div>'
                        f'<div class="metric-value">{len(results)}</div></div>', unsafe_allow_html=True)
        with c2:
            st.markdown(f'<div class="metric-card"><div class="metric-label">Top neglected gene</div>'
                        f'<div class="metric-value" style="color:#fb923c">{results[0]["gene"]}</div></div>',
                        unsafe_allow_html=True)
        with c3:
            top_score = results[0]["score"]
            st.markdown(f'<div class="metric-card"><div class="metric-label">Top neglect score</div>'
                        f'<div class="metric-value" style="color:#fb923c">{top_score}</div></div>',
                        unsafe_allow_html=True)
        with c4:
            avg_papers = int(sum(r["papers"] for r in display) / max(len(display), 1))
            st.markdown(f'<div class="metric-card"><div class="metric-label">Avg papers (top {top_n})</div>'
                        f'<div class="metric-value">{avg_papers}</div></div>', unsafe_allow_html=True)
        with c5:
            cross = sum(1 for r in display if r.get("pathway_count", 1) > 1)
            st.markdown(f'<div class="metric-card"><div class="metric-label">Cross-pathway genes</div>'
                        f'<div class="metric-value">{cross}</div></div>', unsafe_allow_html=True)

        st.markdown("---")

        # Charts
        ch1, ch2 = st.columns(2)
        with ch1:
            st.markdown("#### Neglect Score — Top Genes")
            fig_bar = px.bar(
                df, x="gene", y="score", color="score",
                color_continuous_scale=[[0, "#7c2d00"], [0.5, "#ea580c"], [1, "#fb923c"]],
                labels={"gene": "Gene", "score": "Neglect Score"},
                height=300,
                text=df["score"].apply(lambda s: f"{s:.3f}")
            )
            fig_bar.update_traces(textposition="outside", textfont_color="#fb923c")
            fig_bar.update_layout(
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                font_color="#8b949e", coloraxis_showscale=False,
                margin=dict(l=0, r=0, t=20, b=0),
                xaxis=dict(gridcolor="#21262d"), yaxis=dict(gridcolor="#21262d"),
            )
            st.plotly_chart(fig_bar, use_container_width=True)

        with ch2:
            st.markdown("#### Importance vs. Research Attention")
            fig_scatter = px.scatter(
                df, x="papers", y="centrality",
                size="go_terms", color="score", hover_name="gene",
                color_continuous_scale=[[0, "#7c2d00"], [1, "#fb923c"]],
                labels={"papers": "PubMed Papers", "centrality": "Network Centrality",
                        "go_terms": "GO Terms", "score": "Neglect Score"},
                height=300
            )
            fig_scatter.update_layout(
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                font_color="#8b949e", coloraxis_showscale=False,
                margin=dict(l=0, r=0, t=20, b=0),
                xaxis=dict(gridcolor="#21262d"), yaxis=dict(gridcolor="#21262d"),
            )
            for _, row in df.iterrows():
                fig_scatter.add_annotation(
                    x=row["papers"], y=row["centrality"], text=row["gene"],
                    showarrow=False, font=dict(size=10, color="#c9d1d9"), yshift=10
                )
            st.plotly_chart(fig_scatter, use_container_width=True)

        # Signal breakdown stacked bar
        st.markdown("#### Neglect Score Signal Breakdown")
        st.caption(
            "Stacked bars show the three normalized signals driving each gene's neglect score: "
            "network centrality, GO functional entropy, and inverse literature attention."
        )
        max_c = max((r["centrality"] for r in display), default=1) or 1
        max_e = max((r["entropy"]    for r in display), default=1) or 1
        max_p = max((r["papers"]     for r in display), default=1) or 1
        bdf = pd.DataFrame([{
            "Gene":          r["gene"],
            "Centrality":    round(r["centrality"] / max_c, 3),
            "GO Entropy":    round(r["entropy"]    / max_e, 3),
            "Low Attention": round(1 - r["papers"] / max_p, 3),
        } for r in display])
        fig_bd = go.Figure()
        for signal, color in [("Centrality","#58a6ff"),("GO Entropy","#3fb950"),("Low Attention","#f97316")]:
            fig_bd.add_trace(go.Bar(name=signal, x=bdf["Gene"], y=bdf[signal],
                                    marker_color=color, opacity=0.85))
        fig_bd.update_layout(
            barmode="stack", paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            font_color="#8b949e", legend=dict(font=dict(color="#c9d1d9")),
            margin=dict(l=0, r=0, t=10, b=0), height=300,
            xaxis=dict(gridcolor="#21262d"),
            yaxis=dict(gridcolor="#21262d", title="Normalized signal"),
        )
        st.plotly_chart(fig_bd, use_container_width=True)

        # Cross-pathway heatmap
        if any(r.get("pathway_count", 0) > 1 for r in display):
            st.markdown("#### Cross-Pathway Membership")
            all_pw = sorted(set(p for r in display for p in r.get("pathway_names", [])[:8]))
            matrix = [[1 if p in r.get("pathway_names",[]) else 0 for p in all_pw] for r in display]
            fig_heat = go.Figure(data=go.Heatmap(
                z=matrix, x=[p[:32] for p in all_pw], y=[r["gene"] for r in display],
                colorscale=[[0,"#0d1117"],[1,"#f97316"]], showscale=False,
            ))
            fig_heat.update_layout(
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                font_color="#8b949e", height=max(250, top_n * 28),
                margin=dict(l=0, r=0, t=10, b=60),
                xaxis=dict(tickangle=-35, tickfont=dict(size=10)),
                yaxis=dict(tickfont=dict(size=11, family="IBM Plex Mono")),
            )
            st.plotly_chart(fig_heat, use_container_width=True)

        # Ranking table
        st.markdown("#### Full Gene Ranking")
        table_df = pd.DataFrame([{
            "Gene":          r["gene"],
            "Neglect Score": f"{r['score']:.4f}",
            "PubMed Papers": r["papers"],
            "Centrality":    f"{r['centrality']:.4f}",
            "GO Terms":      r["go_terms"],
            "Entropy":       f"{r['entropy']:.3f}",
            "Trend":         f"{r['trend']:.1f}x",
            "Pathways":      r.get("pathway_count", 1),
        } for r in display])
        st.dataframe(table_df, use_container_width=True, hide_index=True)

        # CSV download
        csv_df = pd.DataFrame(results[:top_n])
        for col in ["pathway_names", "pathway_ids", "years"]:
            if col in csv_df.columns:
                csv_df[col] = csv_df[col].apply(
                    lambda x: "; ".join(str(i) for i in x) if isinstance(x, list) else x
                )
        st.download_button(
            "⬇  Download results CSV",
            csv_df.to_csv(index=False).encode("utf-8"),
            "gene_discovery_results.csv", mime="text/csv"
        )

        # ── AI hypothesis generation ──────────────────────────────────────────
        st.markdown("---")
        st.markdown("## AI Hypothesis Generation")
        st.markdown("Select genes to analyze. Mistral will generate structured research hypotheses.")

        g_col, b_col = st.columns([3, 1])
        with g_col:
            genes_for_ai = st.multiselect(
                "Select genes",
                options=[r["gene"] for r in display],
                default=[display[0]["gene"]] if display else []
            )
        with b_col:
            run_ai = st.button("🤖 Generate Hypotheses", type="primary", use_container_width=True)

        if run_ai and genes_for_ai:
            gene_data_map = {r["gene"]: r for r in results}
            for gene in genes_for_ai:
                if gene in st.session_state.ai_outputs:
                    continue
                with st.spinner(f"Mistral analyzing {gene}..."):
                    data = gene_data_map.get(gene, {"gene": gene})
                    raw  = generate_hypothesis(data, backend="ollama")
                    st.session_state.ai_outputs[gene] = {
                        "raw": raw, "parsed": parse_hypothesis(raw), "data": data
                    }

        # Render AI outputs
        for gene, output in st.session_state.ai_outputs.items():
            parsed = output["parsed"]
            data   = output["data"]
            score  = data.get("score", 0)
            papers = data.get("papers", "?")

            with st.expander(
                f"🧬 {gene}  —  neglect score {score}  |  {papers} papers",
                expanded=True
            ):
                left, right = st.columns([3, 1])

                with left:
                    for section, label, css in [
                        ("biological_importance", "Biological Importance", ""),
                        ("knowledge_gap",          "Knowledge Gap",         "gap"),
                        ("mechanistic_hypothesis", "Mechanistic Hypothesis",""),
                        ("experimental_approach",  "Experimental Approach", "experiment"),
                        ("disease_relevance",       "Disease Relevance",    "disease"),
                        ("drug_targeting",          "Drug Targeting Potential", "drug"),
                    ]:
                        if parsed.get(section):
                            st.markdown(f'<div class="section-header">{label}</div>', unsafe_allow_html=True)
                            st.markdown(f'<div class="hypothesis-section {css}">{parsed[section]}</div>', unsafe_allow_html=True)

                with right:
                    # Percentile rank
                    pct = percentile_rank(score)
                    pct_label = f"top {100 - pct}%" if pct >= 50 else f"bottom {pct + 1}%"

                    # Prominent neglect score + percentile
                    st.markdown(f"""
                    <div class="neglect-card">
                      <div class="nlabel">Neglect Score</div>
                      <div class="nvalue">{score}</div>
                      <div class="nsub">{pct_label} of {len(results)} genes scored</div>
                    </div>""", unsafe_allow_html=True)

                    # Novelty + feasibility with unparsed fallback
                    novelty     = parsed.get("novelty_rating")
                    feasibility = parsed.get("feasibility_rating")
                    parsed_ok   = parsed.get("parsed_ok", False)

                    def rating_card(label, val, gradient):
                        if val is not None:
                            bar = f'<div class="rating-bar-fill" style="width:{val*10}%;{gradient}"></div>'
                            return f"""
                            <div class="metric-card" style="margin-top:6px">
                              <div class="metric-label">{label}</div>
                              <div class="metric-value">{val}<span style="font-size:14px;color:#8b949e">/10</span></div>
                              <div class="rating-bar-bg">{bar}</div>
                            </div>"""
                        else:
                            return f"""
                            <div class="metric-card" style="margin-top:6px;border-color:#374151">
                              <div class="metric-label">{label}</div>
                              <div style="font-size:12px;color:#6b7280;font-family:'IBM Plex Mono',monospace;margin-top:4px">
                                Not parsed — see raw output ↓
                              </div>
                            </div>"""

                    st.markdown(
                        rating_card("Novelty",      novelty,     "background:linear-gradient(90deg,#58a6ff,#3fb950)") +
                        rating_card("Feasibility",  feasibility, "background:linear-gradient(90deg,#22c55e,#4ade80)"),
                        unsafe_allow_html=True
                    )

                    # Publication sparkline (mini)
                    years = data.get("years", [])
                    if years:
                        st.markdown('<div class="section-header">Publication Trend</div>', unsafe_allow_html=True)
                        yc = {}
                        for y in years:
                            yc[y] = yc.get(y, 0) + 1
                        spark_df = pd.DataFrame(sorted(yc.items()), columns=["Year", "Papers"])
                        fig_spark = px.area(spark_df, x="Year", y="Papers",
                            color_discrete_sequence=["#f97316"], height=80)
                        fig_spark.update_layout(
                            paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                            margin=dict(l=0, r=0, t=0, b=0), showlegend=False,
                            xaxis=dict(visible=False), yaxis=dict(visible=False),
                        )
                        fig_spark.update_traces(fillcolor="rgba(249,115,22,0.12)")
                        st.plotly_chart(fig_spark, use_container_width=True, key=f"spark_{gene}")

                    # Signal bars — normalized relative to display set
                    st.markdown('<div class="section-header">Score Signals</div>', unsafe_allow_html=True)
                    centrality = data.get("centrality", 0)
                    entropy    = data.get("entropy",    0)
                    trend      = data.get("trend",      1)
                    c_pct = min(int(centrality / sig_max_c * 100), 100)
                    e_pct = min(int(entropy    / sig_max_e * 100), 100)
                    a_pct = max(0, 100 - min(int(papers / max(sig_max_p, 1) * 100), 100))
                    t_pct = min(int(trend      / sig_max_t * 100), 100)

                    go_note = ' <span style="font-size:10px;color:#6b7280">(no UniProt data)</span>' if entropy == 0 else ""
                    st.markdown(f"""
                    <div class="sig-row">
                      <span class="sig-label">Centrality</span>
                      <div class="sig-bar-bg"><div class="sig-bar-fill" style="width:{c_pct}%;background:#58a6ff"></div></div>
                      <span class="sig-val">{centrality}</span>
                    </div>
                    <div class="sig-row">
                      <span class="sig-label">GO entropy{go_note}</span>
                      <div class="sig-bar-bg"><div class="sig-bar-fill" style="width:{e_pct}%;background:#3fb950"></div></div>
                      <span class="sig-val">{entropy}</span>
                    </div>
                    <div class="sig-row">
                      <span class="sig-label">Low attention</span>
                      <div class="sig-bar-bg"><div class="sig-bar-fill" style="width:{a_pct}%;background:#f97316"></div></div>
                      <span class="sig-val">{papers} pubs</span>
                    </div>
                    <div class="sig-row">
                      <span class="sig-label">Trend</span>
                      <div class="sig-bar-bg"><div class="sig-bar-fill" style="width:{t_pct}%;background:#a855f7"></div></div>
                      <span class="sig-val">{trend}x</span>
                    </div>""", unsafe_allow_html=True)

                    # Pathways — full name in title attr for tooltip
                    pathways = data.get("pathway_names", [])
                    if pathways:
                        st.markdown('<div class="section-header">Pathways</div>', unsafe_allow_html=True)
                        st.markdown(
                            "".join(
                                f'<span class="pathway-tag" title="{p}">{p[:26]}{"…" if len(p)>26 else ""}</span>'
                                for p in pathways[:8]
                            ),
                            unsafe_allow_html=True
                        )

                # Surface raw output prominently if ratings failed to parse
                raw_label = "⚠️ Raw AI output (ratings not parsed)" if not parsed_ok else "Raw AI output"
                with st.expander(raw_label, expanded=not parsed_ok):
                    st.code(output["raw"], language=None)


# ═════════════════════════════════════════════════════════════════════════════
# SINGLE GENE ANALYSIS MODE
# ═════════════════════════════════════════════════════════════════════════════
elif mode == "Single Gene Analysis":
    st.markdown("## Analyze a Single Gene")
    st.markdown("Enter any human gene symbol to compute its neglect score and generate an AI hypothesis.")

    # Ensure session state keys exist for this mode
    if "single_gene_data"   not in st.session_state: st.session_state.single_gene_data   = None
    if "single_gene_ai"     not in st.session_state: st.session_state.single_gene_ai     = None
    if "single_gene_symbol" not in st.session_state: st.session_state.single_gene_symbol = ""

    gene_input = st.text_input("Gene symbol", placeholder="e.g. FOXO3, SIRT1, NRF2, MTOR",
                                value=st.session_state.single_gene_symbol)
    run_single = st.button("▶  Analyze", type="primary")

    # Step 1: score the gene and persist to session state
    if run_single and gene_input:
        gene = gene_input.strip().upper()
        st.session_state.single_gene_symbol = gene
        st.session_state.single_gene_data   = None   # reset on new query
        st.session_state.single_gene_ai     = None
        with st.spinner(f"Scoring {gene} across PubMed, STRING, and UniProt..."):
            result = compute_gene_score(gene, verbose=False)
        if not result:
            st.error(f"Could not score {gene}. Check the gene symbol and try again.")
        else:
            st.session_state.single_gene_data = result

    # Step 2: render scores (persists across reruns because it lives in session state)
    data = st.session_state.single_gene_data
    if data:
        st.markdown("---")

        score_col, signals_col = st.columns([1, 2])
        with score_col:
            st.markdown(f"""
            <div class="neglect-card">
              <div class="nlabel">Neglect Score</div>
              <div class="nvalue">{data["score"]}</div>
              <div class="nsub">higher = more understudied</div>
            </div>""", unsafe_allow_html=True)

        with signals_col:
            m1, m2, m3, m4 = st.columns(4)
            for col, label, val in zip(
                [m1, m2, m3, m4],
                ["PubMed Papers", "Centrality", "GO Terms", "Trend"],
                [data["papers"], data["centrality"], data["go_terms"], f"{data['trend']}x"]
            ):
                with col:
                    st.markdown(f"""
                    <div class="metric-card">
                      <div class="metric-label">{label}</div>
                      <div class="metric-value" style="font-size:20px">{val}</div>
                    </div>""", unsafe_allow_html=True)

        # Publication trend chart
        if data.get("years"):
            st.markdown("---")
            st.markdown("#### Publication Trend")
            year_counts = defaultdict(int)
            for y in data["years"]:
                year_counts[y] += 1
            years_df = pd.DataFrame(sorted(year_counts.items()), columns=["Year", "Papers"])
            fig_trend = px.area(years_df, x="Year", y="Papers",
                color_discrete_sequence=["#f97316"], height=200)
            fig_trend.update_layout(
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                font_color="#8b949e", margin=dict(l=0, r=0, t=10, b=0),
                xaxis=dict(gridcolor="#21262d"), yaxis=dict(gridcolor="#21262d"),
            )
            fig_trend.update_traces(fillcolor="rgba(249,115,22,0.1)")
            st.plotly_chart(fig_trend, use_container_width=True)

        # Step 3: Generate hypothesis button — data is safely in session state
        st.markdown("---")
        st.markdown("#### AI Hypothesis (Mistral)")
        gen_btn = st.button("🤖 Generate Hypothesis", type="primary")

        if gen_btn:
            with st.spinner("Mistral is analyzing..."):
                raw    = generate_hypothesis(data, backend="ollama")
                parsed = parse_hypothesis(raw)
            st.session_state.single_gene_ai = {"raw": raw, "parsed": parsed}

        # Step 4: render AI output (also persists across reruns)
        ai = st.session_state.single_gene_ai
        if ai:
            parsed    = ai["parsed"]
            raw       = ai["raw"]
            parsed_ok = parsed.get("parsed_ok", False)

            left, right = st.columns([3, 1])
            with left:
                for section, label, css in [
                    ("biological_importance", "Biological Importance", ""),
                    ("knowledge_gap",          "Knowledge Gap",         "gap"),
                    ("mechanistic_hypothesis", "Mechanistic Hypothesis",""),
                    ("experimental_approach",  "Experimental Approach", "experiment"),
                    ("disease_relevance",       "Disease Relevance",    "disease"),
                    ("drug_targeting",          "Drug Targeting",       "drug"),
                ]:
                    if parsed.get(section):
                        st.markdown(f'<div class="section-header">{label}</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="hypothesis-section {css}">{parsed[section]}</div>', unsafe_allow_html=True)

            with right:
                novelty     = parsed.get("novelty_rating")
                feasibility = parsed.get("feasibility_rating")

                def _rating_card(label, val, gradient):
                    if val is not None:
                        bar = f'<div class="rating-bar-fill" style="width:{val*10}%;{gradient}"></div>'
                        return f"""
                        <div class="metric-card" style="margin-top:6px">
                          <div class="metric-label">{label}</div>
                          <div class="metric-value">{val}<span style="font-size:14px;color:#8b949e">/10</span></div>
                          <div class="rating-bar-bg">{bar}</div>
                        </div>"""
                    return f"""
                    <div class="metric-card" style="margin-top:6px;border-color:#374151">
                      <div class="metric-label">{label}</div>
                      <div style="font-size:12px;color:#6b7280;font-family:'IBM Plex Mono',monospace;margin-top:4px">
                        Not parsed — see raw output ↓
                      </div>
                    </div>"""

                st.markdown(
                    _rating_card("Novelty",     novelty,     "background:linear-gradient(90deg,#58a6ff,#3fb950)") +
                    _rating_card("Feasibility", feasibility, "background:linear-gradient(90deg,#22c55e,#4ade80)"),
                    unsafe_allow_html=True
                )

                entropy = data.get("entropy", 0)
                if entropy == 0:
                    st.markdown(
                        "<div style='font-size:11px;color:#6b7280;font-family:IBM Plex Mono,monospace;"
                        "margin-top:8px'>⚠ GO entropy = 0: no UniProt data retrieved</div>",
                        unsafe_allow_html=True
                    )

            raw_label = "⚠️ Raw AI output (ratings not parsed)" if not parsed_ok else "Raw AI output"
            with st.expander(raw_label, expanded=not parsed_ok):
                st.code(raw, language=None)