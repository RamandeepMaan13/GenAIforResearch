import requests
import time
from functools import lru_cache

KEGG_BASE = "https://rest.kegg.jp"

@lru_cache(maxsize=256)
def get_kegg_pathway_genes(pathway_id="hsa04115"):
    """
    Retrieve gene symbols from a KEGG pathway.
    Returns an order-preserving deduplicated list.
    """
    try:
        r = requests.get(f"{KEGG_BASE}/get/{pathway_id}", timeout=10)
        time.sleep(0.2)
    except Exception as e:
        print(f"Error fetching pathway {pathway_id}: {e}")
        return []

    if r.status_code != 200:
        print(f"KEGG returned {r.status_code} for {pathway_id}")
        return []

    genes = []
    seen = set()
    reading_genes = False

    for line in r.text.split("\n"):
        if line.startswith("GENE"):
            reading_genes = True
            line = line.replace("GENE", "").strip()
        elif line and not line.startswith(" "):
            reading_genes = False

        if reading_genes:
            parts = line.strip().split()
            if len(parts) >= 2:
                symbol = parts[1].replace(";", "")
                if symbol and symbol not in seen:
                    seen.add(symbol)
                    genes.append(symbol)

    return genes


@lru_cache(maxsize=1)
def get_all_kegg_pathways():
    """Retrieve all human KEGG pathway IDs and names."""
    try:
        r = requests.get(f"{KEGG_BASE}/list/pathway/hsa", timeout=10)
        time.sleep(0.2)
    except Exception as e:
        print(f"Error fetching pathway list: {e}")
        return []

    if r.status_code != 200:
        print(f"KEGG returned {r.status_code}")
        return []

    pathways = []
    for line in r.text.split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            pathway_id = parts[0].replace("path:", "").strip()
            pathway_name = parts[1].strip() if len(parts) > 1 else pathway_id
            pathways.append((pathway_id, pathway_name))

    return pathways


if __name__ == "__main__":
    pathway = input("Enter KEGG pathway ID (example: hsa04115): ")
    genes = get_kegg_pathway_genes(pathway)
    print(f"\nTotal genes retrieved: {len(genes)}")
    print("Sample genes:")
    for g in genes[:10]:
        print(" ", g)