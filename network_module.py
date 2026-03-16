import requests
import time
import networkx as nx
from functools import lru_cache

STRING_BASE = "https://string-db.org/api/json"

@lru_cache(maxsize=512)
def get_string_interactions(gene_name, species=9606, min_score=400):
    """
    Query STRING database for protein interactions.
    min_score: combined score threshold (0-1000). 400 = medium confidence.
    """
    try:
        r = requests.get(f"{STRING_BASE}/network", params={
            "identifiers": gene_name,
            "species": species,
            "required_score": min_score
        }, timeout=15)
        time.sleep(0.3)
    except Exception as e:
        print(f"STRING query failed for {gene_name}: {e}")
        return ()

    if r.status_code != 200:
        return ()

    try:
        data = r.json()
        return tuple((item["preferredName_A"], item["preferredName_B"]) for item in data)
    except Exception:
        return ()


def get_string_interactions_batch(gene_list, species=9606, min_score=400):
    """
    Query STRING for multiple genes in one API call.
    Much faster than one call per gene.
    """
    if not gene_list:
        return []
    try:
        r = requests.get(f"{STRING_BASE}/network", params={
            "identifiers": "%0d".join(gene_list),
            "species": species,
            "required_score": min_score
        }, timeout=20)
        time.sleep(0.3)
    except Exception as e:
        print(f"STRING batch query failed: {e}")
        return []

    if r.status_code != 200:
        return []

    try:
        data = r.json()
        return [(item["preferredName_A"], item["preferredName_B"]) for item in data]
    except Exception:
        return []


def build_graph(interactions):
    """Build undirected graph from interaction pairs."""
    G = nx.Graph()
    for p1, p2 in interactions:
        G.add_edge(p1, p2)
    return G


def compute_centrality(G, gene_name):
    """Compute degree centrality of gene in graph."""
    if gene_name not in G:
        return 0.0
    centrality = nx.degree_centrality(G)
    return centrality.get(gene_name, 0.0)


def get_gene_centrality(gene_name, species=9606):
    """
    Convenience function: fetch interactions and return centrality score.
    No network expansion — fast and reliable.
    """
    interactions = get_string_interactions(gene_name, species)
    if not interactions:
        return 0.0
    G = build_graph(interactions)
    return compute_centrality(G, gene_name)


if __name__ == "__main__":
    gene = input("Enter gene name: ")
    print(f"\nRetrieving interactions for {gene}...")
    score = get_gene_centrality(gene)
    print(f"Network degree centrality: {score:.4f}")