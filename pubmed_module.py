import requests
import time
import re
from xml.etree import ElementTree as ET
from functools import lru_cache

PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
REQUEST_DELAY = 0.34  # ~3 requests/sec without API key

def _get(url, params, retries=3):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=10)
            if r.status_code == 200:
                return r
            time.sleep(1)
        except Exception as e:
            if attempt == retries - 1:
                print(f"Request failed: {e}")
    return None

@lru_cache(maxsize=512)
def get_pubmed_count(gene_name):
    """Get total number of PubMed articles mentioning the gene."""
    r = _get(f"{PUBMED_BASE}/esearch.fcgi", {
        "db": "pubmed",
        "term": f"{gene_name}[Title/Abstract]",
        "retmode": "json"
    })
    time.sleep(REQUEST_DELAY)
    if not r:
        return 0
    try:
        return int(r.json()["esearchresult"]["count"])
    except Exception:
        return 0

@lru_cache(maxsize=512)
def get_publication_years(gene_name):
    """Retrieve publication years for sampled PubMed papers."""
    r = _get(f"{PUBMED_BASE}/esearch.fcgi", {
        "db": "pubmed",
        "term": f"{gene_name}[Title/Abstract]",
        "retmax": 300,
        "retmode": "json"
    })
    time.sleep(REQUEST_DELAY)
    if not r:
        return ()

    try:
        id_list = r.json()["esearchresult"]["idlist"]
        if not id_list:
            return ()
    except Exception:
        return ()

    # Batch fetch in chunks of 100
    years = []
    for i in range(0, len(id_list), 100):
        chunk = id_list[i:i+100]
        fr = _get(f"{PUBMED_BASE}/efetch.fcgi", {
            "db": "pubmed",
            "id": ",".join(chunk),
            "retmode": "xml"
        })
        time.sleep(REQUEST_DELAY)
        if not fr:
            continue
        try:
            root = ET.fromstring(fr.text)
            for article in root.findall(".//PubmedArticle"):
                # Try <Year> first, then <MedlineDate>
                year_el = article.find(".//PubDate/Year")
                if year_el is not None:
                    try:
                        years.append(int(year_el.text))
                        continue
                    except ValueError:
                        pass
                medline_el = article.find(".//PubDate/MedlineDate")
                if medline_el is not None:
                    m = re.search(r'\b(19|20)\d{2}\b', medline_el.text or "")
                    if m:
                        years.append(int(m.group()))
        except ET.ParseError as e:
            print(f"XML parse error for {gene_name}: {e}")

    return tuple(years)

if __name__ == "__main__":
    gene = input("Enter gene name: ")
    count = get_pubmed_count(gene)
    print(f"\nTotal PubMed papers: {count}")
    years = get_publication_years(gene)
    print(f"Papers with year data: {len(years)}")
    if years:
        print(f"Earliest: {min(years)}  Latest: {max(years)}")