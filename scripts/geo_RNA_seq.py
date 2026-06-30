import requests
import pandas as pd
import time
import xml.etree.ElementTree as ET

def ncbi_get(url, params, max_retries=4, timeout=30):
  
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code == 429:
                wait = 2 ** attempt * 2
                print(f"    rate limited, waiting {wait}s...")
                time.sleep(wait)
                continue
            return r
        except (requests.exceptions.ReadTimeout, requests.exceptions.ConnectionError) as e:
            wait = 2 ** attempt * 2
            print(f"    timeout/connection error (attempt {attempt+1}/{max_retries}), waiting {wait}s...")
            time.sleep(wait)
    print(f"    FAILED after {max_retries} attempts: {url}")
    return None

cell_lines = [
    "HL-60","NCI-H929","DND-41","SK-N-SH","Loucy","NB4","MM.1S","MCF7",
    "OCI-LY3","DOHH2","SU-DHL-6","SK-N-MC","LNCaP clone FGC",
    "KMS-11","A549","22Rv1","VCaP","Jurkat"
]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


TREATMENT_KEYWORDS = [
    "treated with", "treatment:", "compound", "drug:", "inhibitor",
    "agonist", "antagonist", "stimulated", "stimulation",
    "transfect", "transduce", "knockdown", "knockout", "sirna", "shrna",
    "crispr", "overexpress", "infect", "irradiat", "hypoxia",
    "starvation", "lps", "cytokine", "dose", "vehicle control",
    "dmso", "ng/ml", "ug/ml", "µg/ml", "nm ", "µm ", " hr ", " hours post"
]

# Subcellular fraction keywords -- exclude anything that isn't whole cell
FRACTION_KEYWORDS = [
    "nuclear fraction", "cytoplasmic fraction", "nucleus",
    "cytoplasm", "mitochondria", "chromatin", "membrane fraction",
    "ribosome", "polysome", "nucleolus", "exosome"
]

def esearch_gds(cell_line):
    """Search GEO DataSets (db=gds) for RNA-seq experiments on this cell line."""
    term = f'"{cell_line}"[All Fields] AND "expression profiling by high throughput sequencing"[DataSet Type]'
    params = {
        "db": "gds",
        "term": term,
        "retmode": "json",
        "retmax": "100"
    }
    r = ncbi_get(f"{EUTILS}/esearch.fcgi", params)
    if r is None or r.status_code != 200:
        return []
    return r.json().get("esearchresult", {}).get("idlist", [])


def esummary_gds(uid_list):
    """Get summaries (GSE accessions, titles) for a list of GDS UIDs."""
    if not uid_list:
        return {}
    params = {
        "db": "gds",
        "id": ",".join(uid_list),
        "retmode": "json"
    }
    r = ncbi_get(f"{EUTILS}/esummary.fcgi", params)
    if r is None or r.status_code != 200:
        return {}
    return r.json().get("result", {})


def get_gse_samples(gse_accession):
    """
    Fetch GSE record from NCBI and parse GSM sample characteristics.
    Returns list of dicts: one per sample with title, characteristics text,
    library_strategy.
    """
    params = {
        "acc": gse_accession,
        "targ": "self",
        "view": "brief",
        "form": "text"
    }
    r = ncbi_get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", params)
    if r is None or r.status_code != 200:
        return None, []

    text = r.text
    library_strategy = None
    for line in text.splitlines():
        if line.lower().startswith("!series_library_strategy") or "library_strategy" in line.lower():
            library_strategy = line.split("=", 1)[-1].strip() if "=" in line else None

    # get sample (GSM) list from the series
    gsm_ids = []
    for line in text.splitlines():
        if line.startswith("!Series_sample_id"):
            gsm_ids.append(line.split("=", 1)[-1].strip())

    return library_strategy, gsm_ids


def get_gsm_characteristics(gsm_accession):
    """Fetch a single GSM record and return its characteristics + treatment text as one string."""
    params = {
        "acc": gsm_accession,
        "targ": "self",
        "view": "brief",
        "form": "text"
    }
    r = ncbi_get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", params)
    if r is None or r.status_code != 200:
        return ""
    return r.text.lower()


def check_untreated_wholecell(sample_text):
    """Returns (untreated: bool, whole_cell: bool) based on keyword scan."""
    treated = any(kw in sample_text for kw in TREATMENT_KEYWORDS)
    fraction = any(kw in sample_text for kw in FRACTION_KEYWORDS)
    return (not treated), (not fraction)



import os

CHECKPOINT_FILE = "geo_rnaseq_filtered_checkpoint.tsv"

rows = []
done_cell_lines = set()

if os.path.exists(CHECKPOINT_FILE):
    prev = pd.read_csv(CHECKPOINT_FILE, sep="\t")
    rows = prev.to_dict("records")
    done_cell_lines = set(prev["cell_line"].unique())
    print(f"Resuming from checkpoint: {len(rows)} rows already collected, "
          f"{len(done_cell_lines)} cell lines done: {sorted(done_cell_lines)}")

for cl in cell_lines:
    if cl in done_cell_lines:
        print(f"\nSkipping {cl} (already done, see checkpoint)")
        continue

    print(f"\nSearching: {cl}")

    uids = esearch_gds(cl)
    print(f"  {len(uids)} GDS/GSE records found")

    if not uids:
        time.sleep(0.4)
        continue

    summaries = esummary_gds(uids)
    time.sleep(0.4)

    for uid in uids:
        info = summaries.get(uid, {})
        gse_acc = info.get("accession", "")
        title = info.get("title", "")

        if not gse_acc.startswith("GSE"):
            print(f"  skip {gse_acc}: not a GSE accession")
            continue

        print(f"  checking {gse_acc}: {title[:60]}")

        library_strategy, gsm_ids = get_gse_samples(gse_acc)
        time.sleep(0.4)

        if library_strategy and "rna-seq" not in library_strategy.lower() \
                and "rnaseq" not in library_strategy.lower():
            print(f"    skip: library_strategy={library_strategy}")
            continue

        if not gsm_ids:
            rows.append({
                "cell_line": cl, "gse_accession": gse_acc, "gsm_accession": "",
                "title": title, "status": "no_samples_found",
                "link": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_acc}"
            })
            continue

        # check each sample (GSM) for treatment / fraction status
        for gsm in gsm_ids[:20]:  
            gsm_text = get_gsm_characteristics(gsm)
            time.sleep(0.34)  
            if "rna-seq" not in gsm_text and "rna seq" not in gsm_text and not library_strategy:
                status = "needs_manual_check"
                untreated, whole_cell = None, None
            else:
                untreated, whole_cell = check_untreated_wholecell(gsm_text)
                if untreated and whole_cell:
                    status = "PASS"
                elif not untreated:
                    status = "treated"
                elif not whole_cell:
                    status = "subcellular_fraction"
                else:
                    status = "needs_manual_check"

            rows.append({
                "cell_line": cl,
                "gse_accession": gse_acc,
                "gsm_accession": gsm,
                "title": title,
                "status": status,
                "link": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
            })
            print(f"    {gsm}: {status}")

    # checkpoint after each cell line finishes
    pd.DataFrame(rows).to_csv(CHECKPOINT_FILE, sep="\t", index=False)
    print(f"  checkpoint saved: {len(rows)} total rows so far")

    time.sleep(0.5)


df = pd.DataFrame(rows)
df.to_csv("geo_rnaseq_filtered.tsv", sep="\t", index=False)

print("\n=== SUMMARY ===")
print(df["status"].value_counts().to_string())
print()
print("=== PASSED (untreated + whole cell) ===")
passed = df[df["status"] == "PASS"]
print(passed[["cell_line", "gse_accession", "gsm_accession", "title"]].to_string(index=False))
print(f"\nTotal passed: {len(passed)}")
