import math
from pubmed_module import get_pubmed_count, get_publication_years
from network_module import get_gene_centrality
from go_module import get_go_terms, compute_functional_entropy

def compute_attention_trend(years):
    """
    Compare papers in last 5 years vs previous 5 years.
    Returns ratio > 1 if research is growing, < 1 if declining.
    """
    if not years:
        return 1.0
    current_year = max(years)
    recent = [y for y in years if y >= current_year - 5]
    previous = [y for y in years if current_year - 10 <= y < current_year - 5]
    if len(previous) == 0:
        return 2.0 if recent else 1.0
    return len(recent) / len(previous)


def compute_gene_score(gene, verbose=False):
    """
    Compute neglect score for a gene.
    
    Returns a dict with all signals and the final score,
    or None if scoring fails entirely.
    
    Score formula:
        importance = normalize(centrality) + normalize(entropy) + trend_bonus
        neglect_score = importance / log(papers + 2)
    
    Higher score = more biologically important but understudied.
    """
    if verbose:
        print(f"\n{'='*40}\nAnalyzing: {gene}\n{'='*40}")

    try:
        # 1. Literature attention
        papers = get_pubmed_count(gene)
        if verbose:
            print(f"PubMed papers: {papers}")

        # 2. Attention trend
        years = get_publication_years(gene)
        trend = compute_attention_trend(years)
        if verbose:
            print(f"Attention trend (recent/previous 5yr): {trend:.2f}")

        # 3. Network centrality (no expand — fast)
        centrality = get_gene_centrality(gene)
        if verbose:
            print(f"Network centrality: {centrality:.4f}")

        # 4. Functional diversity
        go_terms = get_go_terms(gene)
        entropy = compute_functional_entropy(go_terms)
        if verbose:
            print(f"Functional diversity (GO entropy): {entropy:.4f}")
            print(f"GO terms: {len(go_terms)}")

        # Normalize centrality to match entropy scale
        # centrality: 0-1, entropy: 0 ~ log(3)*log(n+1) ~ up to ~5
        # Bring centrality into comparable range
        centrality_scaled = centrality * math.log(len(go_terms) + 2)

        # Trend bonus: growing research area gets a small boost (0 to 0.3)
        trend_bonus = min(0.3, math.log(trend + 1) * 0.2)

        importance = centrality_scaled + entropy + trend_bonus

        # Penalize by research attention (log-scaled)
        neglect_score = importance / math.log(papers + 2)

        if verbose:
            print(f"\nImportance score: {importance:.4f}")
            print(f"Final neglect score: {neglect_score:.4f}")

        return {
            "gene": gene,
            "score": round(neglect_score, 4),
            "papers": papers,
            "centrality": round(centrality, 4),
            "go_terms": len(go_terms),
            "entropy": round(entropy, 4),
            "trend": round(trend, 2),
            "years": list(years) if years else []
        }

    except Exception as e:
        print(f"Scoring failed for {gene}: {e}")
        return None


if __name__ == "__main__":
    gene = input("Enter gene name: ")
    result = compute_gene_score(gene, verbose=True)
    if result:
        print(f"\nResult: {result}")