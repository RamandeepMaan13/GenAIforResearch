import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from pathway_gene_fetcher import get_kegg_pathway_genes, get_all_kegg_pathways
from gene_score_engine import compute_gene_score

def scan_pathways(
    max_pathways=None,
    max_workers=6,
    progress_callback=None,
    pathway_filter=None
):
    """
    Scan KEGG human pathways and score all unique genes.

    Args:
        max_pathways: limit number of pathways (None = all)
        max_workers: parallel threads for gene scoring
        progress_callback: optional fn(current, total, message) for UI updates
        pathway_filter: optional fn(pathway_name) -> bool to pre-filter pathways
                        before scanning. Filtering happens before any gene API calls
                        so irrelevant pathways are completely skipped.

    Returns:
        List of dicts: {gene, score, papers, centrality, go_terms, entropy, trend, pathways}
    """
    all_pathways = get_all_kegg_pathways()

    # Apply theme filter BEFORE slicing — so max_pathways applies to the matched set
    if pathway_filter:
        pathways = [(pid, pname) for pid, pname in all_pathways if pathway_filter(pname)]
        if progress_callback:
            progress_callback(0, len(pathways),
                f"Theme filter: {len(pathways)} / {len(all_pathways)} pathways matched")
    else:
        pathways = list(all_pathways)

    if max_pathways:
        pathways = pathways[:max_pathways]

    total_pathways = len(pathways)
    if not total_pathways:
        if progress_callback:
            progress_callback(0, 0, "No pathways matched the theme filter — try a different theme")
        return []

    if progress_callback:
        progress_callback(0, total_pathways, f"Scanning {total_pathways} pathways...")

    # Step 1: collect all genes and their pathway memberships
    gene_pathway_map = defaultdict(list)  # gene -> [(pathway_id, pathway_name)]

    for i, (pathway_id, pathway_name) in enumerate(pathways):
        genes = get_kegg_pathway_genes(pathway_id)
        for gene in genes:
            gene_pathway_map[gene].append((pathway_id, pathway_name))
        if progress_callback:
            progress_callback(i + 1, total_pathways,
                f"Scanning pathway {i+1}/{total_pathways}: {pathway_name} ({len(genes)} genes)")

    unique_genes = list(gene_pathway_map.keys())
    total_genes = len(unique_genes)

    if progress_callback:
        progress_callback(0, total_genes,
            f"Found {total_genes} unique genes — scoring now...")

    # Step 2: score each unique gene in parallel
    scored = {}
    completed = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {executor.submit(compute_gene_score, gene): gene for gene in unique_genes}
        for future in as_completed(future_map):
            gene = future_map[future]
            completed += 1
            try:
                result = future.result()
                if result is not None:
                    scored[gene] = result
            except Exception as e:
                print(f"Scoring failed for {gene}: {e}")
            if progress_callback:
                progress_callback(completed, total_genes,
                    f"Scored {completed}/{total_genes} genes — last: {gene}")

    # Step 3: merge pathway info into results
    results = []
    for gene, data in scored.items():
        pathway_entries = gene_pathway_map[gene]
        pathway_ids = [p for p, _ in pathway_entries]
        pathway_names = [n for _, n in pathway_entries]
        cross_pathway_boost = 0.1 * (len(pathway_ids) - 1)
        final_score = round(data["score"] + cross_pathway_boost, 4)

        results.append({
            **data,
            "score": final_score,
            "base_score": data["score"],
            "pathway_ids": pathway_ids,
            "pathway_names": pathway_names,
            "pathway_count": len(pathway_ids)
        })

    results.sort(key=lambda x: x["score"], reverse=True)
    return results


if __name__ == "__main__":
    def cli_progress(current, total, message):
        print(f"  [{current}/{total}] {message}")

    results = scan_pathways(max_pathways=10, progress_callback=cli_progress)
    print(f"\n{'='*50}")
    print("TOP NEGLECTED GENES")
    print(f"{'='*50}\n")
    for r in results[:15]:
        print(f"{r['gene']:10s} score={r['score']:.4f}  papers={r['papers']}  pathways={r['pathway_count']}")