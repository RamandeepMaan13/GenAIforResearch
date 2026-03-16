import requests
import csv
from collections import defaultdict
from global_pathway_scan import scan_pathways


# ---------------- AI ENGINE ---------------- #

def run_ai_analysis(gene, pathways, score):

    pathway_text = ", ".join(pathways)

    prompt = f"""
You are a senior molecular biologist.

Gene: {gene}
Associated pathways: {pathway_text}
Neglect score: {score}

Provide structured output:

1. Biological importance (2 lines max)
2. Key knowledge gap
3. Mechanistic hypothesis
4. Experimental approach
5. Disease relevance (if any)
6. Drug targeting potential (high/medium/low + reason)
7. Novelty rating (1–10)
8. Feasibility rating (1–10)

Be concise and scientific.
"""

    response = requests.post(
        "http://localhost:11434/api/generate",
        json={
            "model": "mistral",
            "prompt": prompt,
            "stream": False
        }
    )

    return response.json()["response"]


# ---------------- CROSS-PATHWAY DETECTOR ---------------- #

def detect_cross_pathway_genes(ranked):

    gene_map = defaultdict(list)

    for gene, pathway, score in ranked:
        gene_map[gene].append((pathway, score))

    results = []

    for gene, entries in gene_map.items():

        pathways = [p for p, s in entries]
        avg_score = sum(s for p, s in entries) / len(entries)

        novelty_boost = len(pathways)

        final_score = avg_score + (0.1 * novelty_boost)

        results.append((gene, pathways, final_score))

    return sorted(results, key=lambda x: x[2], reverse=True)


# ---------------- EXPORT SYSTEM ---------------- #

def export_results(rows):

    with open("final_translational_discovery_results.csv", "w", newline="") as f:

        writer = csv.writer(f)

        writer.writerow(["Gene", "Pathways", "Score", "AI Output"])

        for row in rows:
            writer.writerow(row)


# ---------------- MAIN SYSTEM ---------------- #

if __name__ == "__main__":

    print("\nScanning biology...\n")

    ranked = scan_pathways()

    print("\nDetecting cross-pathway discoveries...\n")

    cross_ranked = detect_cross_pathway_genes(ranked)

    print("\nRunning AI reasoning...\n")

    results = []

    for gene, pathways, score in cross_ranked[:5]:

        print(f"\nGene: {gene} | Score: {score}")

        ai_output = run_ai_analysis(gene, pathways, score)

        print("\n----------------")
        print(ai_output)
        print("----------------\n")

        results.append([gene, ", ".join(pathways), score, ai_output])

    export_results(results)

    print("\nSaved to final_translational_discovery_results.csv\n")