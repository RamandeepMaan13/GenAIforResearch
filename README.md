# 🧬 Gene Discovery Assistant

An AI-powered bioinformatics tool that identifies **neglected genes** — genes with high biological importance that are significantly understudied relative to their research potential. It combines multi-source biological signals with a local LLM (Mistral via Ollama) to generate structured research hypotheses.

---

## What It Does

The tool scores every gene in human KEGG pathways using three independent signals:

- **Network centrality** — how connected the gene is in the STRING protein interaction network
- **Functional diversity** — how broadly the gene is annotated across GO categories (via UniProt)
- **Literature attention** — how many PubMed papers mention the gene

Genes that score high on importance but low on attention rise to the top as neglected candidates. Mistral then generates structured hypotheses for each: biological importance, knowledge gaps, mechanistic hypotheses, experimental approaches, disease relevance, and drug targeting potential.

---

## Project Structure

```
├── app.py                      # Streamlit frontend
├── gene_score_engine.py        # Core neglect score computation
├── global_pathway_scan.py      # Scans all human KEGG pathways in parallel
├── pathway_gene_fetcher.py     # Fetches gene lists from KEGG API
├── pubmed_module.py            # PubMed article counts and publication years
├── network_module.py           # STRING protein interaction network + centrality
├── go_module.py                # UniProt GO terms and functional entropy
├── ai_engine.py                # Unified Ollama/Mistral hypothesis generator
├── pathway_discovery_engine.py # Single pathway analysis
├── gene_ranker.py              # Rank a custom gene list
├── discovery_assistant.py      # Orchestrator + CSV export
├── figure_generator.py         # Matplotlib visualization export
└── requirements.txt
```

---

## Neglect Score Formula

```
importance  = normalize(centrality) + normalize(GO_entropy) + trend_bonus
neglect_score = importance / log(pubmed_papers + 2)
```

Genes appearing in multiple pathways receive a cross-pathway boost. Higher score = more biologically important but understudied.

---

## Prerequisites

- Python 3.10 or higher
- [Ollama](https://ollama.com) installed and running locally

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/RamandeepMaan13/GenAIforResearch.git
cd GenAIforResearch

# 2. (Recommended) Create a virtual environment
python3 -m venv venv
source venv/bin/activate        # Mac / Linux
venv\Scripts\activate           # Windows

# 3. Install Python dependencies
pip install -r requirements.txt

# 4. Pull the Mistral model via Ollama
ollama pull mistral

# 5. Start Ollama (keep this running in a separate terminal)
ollama serve
```

---

## Running the App

```bash
streamlit run app.py
```

Then open [http://localhost:8501](http://localhost:8501) in your browser.

---

## Usage

### Discovery Scan
1. Select a **Biological Theme** from the sidebar (Cancer, Metabolism, Aging & Longevity, etc.) or leave as "All pathways"
2. Set the number of pathways to scan and parallel workers
3. Click **▶ Run Discovery Scan**
4. Watch genes being scored in real time in the live log
5. Explore the ranked results, charts, and cross-pathway heatmap
6. Select genes and click **🤖 Generate Hypotheses** to get Mistral's analysis

### Single Gene Analysis
1. Switch to **Single Gene Analysis** mode in the sidebar
2. Enter any human gene symbol (e.g. `FOXO3`, `SIRT1`, `NRF2`)
3. Click **▶ Analyze** to compute the neglect score
4. Click **🤖 Generate Hypothesis** to get the AI analysis

### Custom Gene List
1. In Discovery Scan mode, set Gene Source to **Custom gene list**
2. Enter comma-separated gene symbols
3. Run the scan

---

## External APIs Used

| API | Purpose | Rate limit |
|-----|---------|------------|
| [KEGG REST API](https://rest.kegg.jp) | Pathway gene lists | Free, no key needed |
| [PubMed E-utilities](https://eutils.ncbi.nlm.nih.gov) | Literature counts | 3 req/sec without API key |
| [STRING DB](https://string-db.org) | Protein interactions | Free, no key needed |
| [UniProt REST](https://rest.uniprot.org) | GO term annotations | Free, no key needed |
| [Ollama](https://ollama.com) | Local LLM inference | Runs on your machine |

No API keys are required. Everything runs locally or against free public APIs.

---

## Biological Themes Available

| Theme | Example KEGG pathways matched |
|-------|-------------------------------|
| Cancer | Colorectal cancer, Glioma, Melanoma |
| Metabolism | Glycolysis, Fatty acid degradation, TCA cycle |
| Aging & Longevity | FoxO signaling, mTOR signaling, Cellular senescence |
| Immune System | T cell receptor signaling, NF-kappa B, Cytokine-cytokine interaction |
| Neurodegeneration | Alzheimer disease, Parkinson disease, ALS |
| Cardiovascular | Dilated cardiomyopathy, Hypertrophic cardiomyopathy |
| Cell Cycle & Apoptosis | Cell cycle, Apoptosis, p53 signaling, Ferroptosis |
| Diabetes & Obesity | Insulin signaling, Adipogenesis, PPAR signaling |
| Inflammation | JAK-STAT, IL-17, TNF, NOD-like receptor |
| Stem Cells & Development | Wnt, Notch, Hedgehog, TGF-beta |
| Infectious Disease | Tuberculosis, HIV, Influenza, Hepatitis |

---

## Performance Tips

- Start with **10–20 pathways** for a quick test run (~2–3 minutes)
- Use **5–6 parallel workers** for a good balance of speed and API stability
- Results are cached in session — switching tabs or adjusting visualizations does not re-trigger the scan
- For a full scan of all 330 human KEGG pathways, allow 30–60 minutes depending on your connection

---

## Output

Results can be downloaded as a CSV containing:

- Gene symbol
- Neglect score
- PubMed paper count
- Network centrality
- GO terms count
- Functional entropy
- Publication trend (recent vs previous 5 years)
- Associated pathway names

---

## Built With

- [Streamlit](https://streamlit.io) — frontend
- [Plotly](https://plotly.com) — interactive charts
- [NetworkX](https://networkx.org) — graph analysis
- [Ollama](https://ollama.com) + [Mistral](https://mistral.ai) — local LLM inference
- KEGG, PubMed, STRING, UniProt — biological databases

---

## Author

Ramandeep Kaur — [GitHub](https://github.com/RamandeepMaan13)
