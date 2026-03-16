import requests
import math
import time
from functools import lru_cache

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

GO_CATEGORIES = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}

@lru_cache(maxsize=512)
def get_go_terms(gene):
    """
    Fetch GO terms from UniProt (Swiss-Prot reviewed entries only).
    Returns a tuple of (go_id, category) pairs.
    """
    try:
        r = requests.get(f"{UNIPROT_BASE}/search", params={
            "query": f"gene:{gene} AND organism_id:9606 AND reviewed:true",
            "format": "json",
            "fields": "xref_go",
            "size": 1
        }, timeout=10)
        time.sleep(0.2)
    except Exception as e:
        print(f"UniProt query failed for {gene}: {e}")
        return ()

    if r.status_code != 200:
        return ()

    go_terms = []
    try:
        data = r.json()
        results = data.get("results", [])
        if not results:
            return ()
        entry = results[0]
        for ref in entry.get("uniProtKBCrossReferences", []):
            if ref.get("database") == "GO":
                go_id = ref.get("id", "")
                # Extract category from GO term properties
                category = "unknown"
                for prop in ref.get("properties", []):
                    if prop.get("key") == "GoTerm":
                        val = prop.get("value", "")
                        if val.startswith("F:"):
                            category = "molecular_function"
                        elif val.startswith("P:"):
                            category = "biological_process"
                        elif val.startswith("C:"):
                            category = "cellular_component"
                if go_id:
                    go_terms.append((go_id, category))
    except Exception as e:
        print(f"GO term parsing failed for {gene}: {e}")

    return tuple(go_terms)


def compute_functional_entropy(go_terms):
    """
    Compute Shannon entropy over GO categories (MF/BP/CC distribution).
    A gene active in all three categories scores higher than one concentrated in one.
    Returns value in [0, log(3)] ~ [0, 1.099].
    """
    if not go_terms:
        return 0.0

    counts = {"molecular_function": 0, "biological_process": 0, "cellular_component": 0}
    for _, cat in go_terms:
        if cat in counts:
            counts[cat] += 1

    total = sum(counts.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log(p)

    # Also factor in sheer number of terms (log-scaled breadth)
    breadth = math.log(len(go_terms) + 1)
    return entropy * breadth


def get_go_summary(gene):
    """Return a dict with GO counts by category for display."""
    terms = get_go_terms(gene)
    summary = {"molecular_function": 0, "biological_process": 0, "cellular_component": 0, "total": len(terms)}
    for _, cat in terms:
        if cat in summary:
            summary[cat] += 1
    return summary


if __name__ == "__main__":
    gene = input("Enter gene name: ")
    terms = get_go_terms(gene)
    print(f"GO terms retrieved: {len(terms)}")
    entropy = compute_functional_entropy(terms)
    print(f"Functional diversity score: {entropy:.4f}")
    summary = get_go_summary(gene)
    print(f"Breakdown: MF={summary['molecular_function']} BP={summary['biological_process']} CC={summary['cellular_component']}")