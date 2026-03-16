from pathway_gene_fetcher import get_kegg_pathway_genes
from gene_score_engine import compute_gene_score


def analyze_pathway(pathway_id):

    print("\nFetching genes from pathway:", pathway_id)

    genes = get_kegg_pathway_genes(pathway_id)

    print("Number of genes retrieved:", len(genes))

    results = []

    for gene in genes:

        try:
            score = compute_gene_score(gene)

            if score is not None:
                results.append((gene, score))

        except Exception as e:
            print("Skipping gene:", gene)

    ranked = sorted(results, key=lambda x: x[1], reverse=True)

    return ranked


if __name__ == "__main__":

    pathway = input("Enter KEGG pathway ID (example: hsa04115): ")

    ranking = analyze_pathway(pathway)

    print("\n==============================")
    print("TOP NEGLECTED GENES IN PATHWAY")
    print("==============================\n")

    for gene, score in ranking[:10]:
        print(gene, ":", score)
