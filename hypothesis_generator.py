import requests


def generate_hypothesis(gene):

    prompt = f"""


Gene: {gene}

Generate:

1. One novel mechanistic hypothesis
2. One research gap
3. One experiment to test it

Be concise and realistic.
"""

    response = requests.post(
        "http://localhost:11434/api/generate",
        json={
            "model": "mistral",
            "prompt": prompt,
            "stream": False
        }
    )

    result = response.json()

    print("\n==============================")
    print(result["response"])
    print("==============================\n")


if __name__ == "__main__":

    gene = input("Enter gene: ")
    generate_hypothesis(gene)