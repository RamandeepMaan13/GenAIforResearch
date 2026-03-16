"""
ai_engine.py — Unified AI backend for gene hypothesis generation.
Supports: OpenAI (GPT-4.1-mini), Ollama (local Mistral), or mock mode.
"""
import os
import requests

# ── Backend selection ──────────────────────────────────────────────────────────
BACKEND = os.environ.get("AI_BACKEND", "ollama")   # "openai" | "ollama" | "mock"
OLLAMA_URL = os.environ.get("OLLAMA_URL", "http://localhost:11434")
OLLAMA_MODEL = os.environ.get("OLLAMA_MODEL", "mistral")
OPENAI_MODEL = os.environ.get("OPENAI_MODEL", "gpt-4.1-mini")


def _build_prompt(gene_data: dict) -> str:
    gene = gene_data.get("gene", "unknown")
    papers = gene_data.get("papers", "N/A")
    score = gene_data.get("score", "N/A")
    centrality = gene_data.get("centrality", "N/A")
    go_terms = gene_data.get("go_terms", "N/A")
    pathways = gene_data.get("pathway_names", [])
    pathway_text = ", ".join(pathways[:5]) if pathways else "unknown"

    return f"""You are a senior molecular biologist reviewing computational gene prioritization results.

Gene: {gene}
Neglect score: {score} (higher = more understudied relative to importance)
PubMed papers: {papers}
Network centrality: {centrality}
GO terms: {go_terms}
Associated pathways: {pathway_text}

Provide a structured analysis with EXACTLY these 8 sections:

1. BIOLOGICAL IMPORTANCE
Why this gene matters biologically (2 sentences max).

2. KNOWLEDGE GAP
The most critical thing we do not yet understand about this gene.

3. MECHANISTIC HYPOTHESIS
One specific, testable mechanistic hypothesis about this gene's function.

4. EXPERIMENTAL APPROACH
The single best experiment to test this hypothesis (be specific: method, model system).

5. DISEASE RELEVANCE
Diseases or conditions this gene may be implicated in (or "None identified").

6. DRUG TARGETING POTENTIAL
High / Medium / Low — with a one-sentence reason.

7. NOVELTY RATING
X/10 — with a one-sentence justification.

8. FEASIBILITY RATING
X/10 — with a one-sentence justification.

Be concise, specific, and scientific. Avoid generic statements."""


def _call_ollama(prompt: str) -> str:
    try:
        r = requests.post(f"{OLLAMA_URL}/api/generate", json={
            "model": OLLAMA_MODEL,
            "prompt": prompt,
            "stream": False
        }, timeout=120)
        r.raise_for_status()
        return r.json().get("response", "").strip()
    except requests.exceptions.ConnectionError:
        return "ERROR: Ollama is not running. Start it with `ollama serve` or switch to OpenAI backend."
    except Exception as e:
        return f"ERROR: Ollama request failed — {e}"


def _call_openai(prompt: str) -> str:
    try:
        from openai import OpenAI
        client = OpenAI()
        response = client.chat.completions.create(
            model=OPENAI_MODEL,
            messages=[{"role": "user", "content": prompt}],
            max_tokens=800,
            temperature=0.7
        )
        return response.choices[0].message.content.strip()
    except ImportError:
        return "ERROR: openai package not installed. Run: pip install openai"
    except Exception as e:
        return f"ERROR: OpenAI request failed — {e}"


def _call_mock(gene_data: dict) -> str:
    gene = gene_data.get("gene", "UNKNOWN")
    return f"""1. BIOLOGICAL IMPORTANCE
{gene} is a computationally prioritized gene with high network centrality relative to its publication count. It participates in multiple signaling cascades.

2. KNOWLEDGE GAP
The precise regulatory mechanisms governing {gene} expression under metabolic stress remain poorly characterized.

3. MECHANISTIC HYPOTHESIS
{gene} may act as a molecular switch linking nutrient sensing to transcriptional reprogramming via post-translational modification of downstream effectors.

4. EXPERIMENTAL APPROACH
CRISPR knockout of {gene} in HEK293T cells followed by phosphoproteomic profiling under glucose deprivation to map downstream signaling changes.

5. DISEASE RELEVANCE
Potential relevance to metabolic syndrome and cancer cachexia based on pathway co-membership.

6. DRUG TARGETING POTENTIAL
Medium — the protein contains a predicted allosteric pocket but lacks validated small-molecule binders.

7. NOVELTY RATING
7/10 — high network importance with relatively few mechanistic studies.

8. FEASIBILITY RATING
8/10 — standard molecular tools available; CRISPR lines readily generated."""


def generate_hypothesis(gene_data: dict, backend: str = None) -> str:
    """
    Generate a structured hypothesis for a gene.
    
    Args:
        gene_data: dict from compute_gene_score()
        backend: override BACKEND env var ("openai", "ollama", "mock")
    
    Returns:
        Formatted text with 8 structured sections.
    """
    active_backend = backend or BACKEND
    prompt = _build_prompt(gene_data)

    if active_backend == "openai":
        return _call_openai(prompt)
    elif active_backend == "mock":
        return _call_mock(gene_data)
    else:
        return _call_ollama(prompt)


def parse_hypothesis(text: str) -> dict:
    """
    Robustly parse structured AI output into a dict.
    Handles many Mistral output variants:
      - "7/10", "7 / 10", "7 out of 10", "**7**/10", "Score: 7", "Rating: 7"
      - Section headers with/without numbers, bold markdown, colons
      - Numbered "1." or "1)" prefixes
    """
    import re

    sections = {
        "biological_importance": "",
        "knowledge_gap":          "",
        "mechanistic_hypothesis": "",
        "experimental_approach":  "",
        "disease_relevance":      "",
        "drug_targeting":         "",
        "novelty_rating":         None,
        "feasibility_rating":     None,
        "parsed_ok":              False,
        "raw":                    text,
    }

    # Map of header fragment → section key (order matters — checked top-to-bottom)
    HEADER_MAP = [
        (["BIOLOGICAL IMPORTANCE", "BIOLOGICAL ROLE", "WHY THIS GENE"],  "biological_importance"),
        (["KNOWLEDGE GAP", "RESEARCH GAP", "KEY GAP"],                   "knowledge_gap"),
        (["MECHANISTIC HYPOTHESIS", "HYPOTHESIS"],                        "mechanistic_hypothesis"),
        (["EXPERIMENTAL APPROACH", "EXPERIMENT"],                         "experimental_approach"),
        (["DISEASE RELEVANCE", "DISEASE"],                                "disease_relevance"),
        (["DRUG TARGETING", "DRUG POTENTIAL", "DRUGGABILITY"],            "drug_targeting"),
        (["NOVELTY RATING", "NOVELTY SCORE", "NOVELTY"],                  "novelty_rating"),
        (["FEASIBILITY RATING", "FEASIBILITY SCORE", "FEASIBILITY"],      "feasibility_rating"),
    ]

    def detect_section(line: str):
        """Return section key if this line is a section header, else None."""
        # Strip markdown bold/italic, leading numbers+punctuation, colons
        clean = re.sub(r'\*+', '', line).strip()
        clean = re.sub(r'^[\d]+[.)]\s*', '', clean).strip()
        upper = clean.upper()
        for fragments, key in HEADER_MAP:
            if any(frag in upper for frag in fragments):
                # Make sure it's short enough to be a header (not a sentence containing the word)
                if len(clean) < 80:
                    return key
        return None

    def flush_buffer(key, buf):
        """Convert buffer lines to section content."""
        content = "\n".join(buf).strip()
        # Remove leading colon or dash
        content = re.sub(r'^[\s:–—-]+', '', content).strip()
        if not content:
            return

        if key in ("novelty_rating", "feasibility_rating"):
            # Try many formats:
            # "7/10", "7 / 10", "**7**/10", "7 out of 10", "Rating: 7", "Score: 7", standalone "7"
            rating = None
            patterns = [
                r'\b([1-9]|10)\s*/\s*10',           # 7/10
                r'\b([1-9]|10)\s+out\s+of\s+10',    # 7 out of 10
                r'rating[:\s]+\*{0,2}([1-9]|10)\*{0,2}',  # Rating: 7
                r'score[:\s]+\*{0,2}([1-9]|10)\*{0,2}',   # Score: 7
                r'^\*{0,2}([1-9]|10)\*{0,2}[\s/–—-]',     # starts with the number
                r'\*{0,2}([1-9]|10)\*{0,2}\s*/\s*10',      # **7**/10
            ]
            for pat in patterns:
                m = re.search(pat, content, re.IGNORECASE)
                if m:
                    rating = int(m.group(1))
                    break
            # Last resort: first standalone digit 1-10 on its own line
            if rating is None:
                for l in content.split("\n"):
                    m = re.match(r'^\s*\*{0,2}(10|[1-9])\*{0,2}\s*$', l.strip())
                    if m:
                        rating = int(m.group(1))
                        break
            sections[key] = rating
        else:
            sections[key] = content

    current_key = None
    buffer = []

    for line in text.split("\n"):
        key = detect_section(line)
        if key:
            # Flush previous buffer
            if current_key:
                flush_buffer(current_key, buffer)
            current_key = key
            buffer = []
            # If the header line itself contains content after a colon, keep it
            after_colon = re.split(r':\s+', line, maxsplit=1)
            if len(after_colon) > 1 and len(after_colon[1].strip()) > 4:
                buffer.append(after_colon[1].strip())
        elif current_key:
            buffer.append(line)

    # Flush final buffer
    if current_key:
        flush_buffer(current_key, buffer)

    # Mark as successfully parsed if at least 3 sections filled
    filled = sum(1 for k, v in sections.items()
                 if k not in ("raw", "parsed_ok") and v)
    sections["parsed_ok"] = filled >= 3

    return sections


if __name__ == "__main__":
    import json
    gene_data = {
        "gene": input("Enter gene: "),
        "score": 0.85,
        "papers": 42,
        "centrality": 0.31,
        "go_terms": 18,
        "pathway_names": ["p53 signaling", "apoptosis"]
    }
    result = generate_hypothesis(gene_data)
    print("\n" + "="*50)
    print(result)
    print("="*50)
    parsed = parse_hypothesis(result)
    print("\nParsed novelty rating:", parsed["novelty_rating"])
    print("Parsed feasibility rating:", parsed["feasibility_rating"])