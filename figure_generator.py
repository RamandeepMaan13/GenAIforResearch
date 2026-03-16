import pandas as pd
import matplotlib.pyplot as plt


# -------- LOAD DATA -------- #

df = pd.read_csv("final_translational_discovery_results.csv")


# -------- FIGURE 1: TOP GENES -------- #

plt.figure()

plt.bar(df["Gene"], df["Score"])

plt.title("Top Neglected Genes Across Pathways")
plt.ylabel("Discovery Score")
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig("figure_top_genes.png")


# -------- FIGURE 2: PATHWAY COUNTS -------- #

pathway_counts = df["Pathways"].apply(lambda x: len(str(x).split(",")))

plt.figure()

plt.bar(df["Gene"], pathway_counts)

plt.title("Cross-Pathway Involvement")
plt.ylabel("Number of Pathways")
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig("figure_cross_pathway.png")


# -------- FIGURE 3: NOVELTY VS FEASIBILITY -------- #

# Extract ratings from AI output
def extract_value(text, keyword):

    try:
        for line in text.split("\n"):
            if keyword in line:
                return float("".join([c for c in line if c.isdigit()]))
    except:
        return None


df["Novelty"] = df["AI Output"].apply(lambda x: extract_value(x, "Novelty"))
df["Feasibility"] = df["AI Output"].apply(lambda x: extract_value(x, "Feasibility"))

plt.figure()

plt.scatter(df["Novelty"], df["Feasibility"])

for i, gene in enumerate(df["Gene"]):
    plt.text(df["Novelty"][i], df["Feasibility"][i], gene)

plt.xlabel("Novelty Score")
plt.ylabel("Feasibility Score")
plt.title("Novelty vs Feasibility")

plt.tight_layout()
plt.savefig("figure_novelty_feasibility.png")

print("\nFigures generated successfully!\n")
