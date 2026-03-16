
from gene_score_engine import compute_gene_score


def rank_genes(gene_list):

    results = []

    for gene in gene_list:

        try:

            score = compute_gene_score(gene)

            results.append((gene, score))

        except Exception as e:

            print("Error analyzing", gene)

    ranked = sorted(results, key=lambda x: x[1], reverse=True)

    return ranked


if __name__ == "__main__":

    genes = [
        "TP53",
        "BRCA1",
        "AKT1",
        "FOXO3",
        "SIRT1",
        "NRF2",
        "MTOR",
        "HIF1A",
        "AMPK",
        "MYC"
    ]

    ranking = rank_genes(genes)

    print("\n==============================")
    print("GENE NEGLECT RANKING")
    print("==============================\n")

    for gene, score in ranking:

        print(gene, ":", score)