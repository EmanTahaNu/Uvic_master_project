

import time
import requests
import pandas as pd

# cell lines NOT found on ENCODE RNA-seq - update this list as your ENCODE
missing_cell_lines = [
    "HL-60", "NCI-H929", "DND-41", "SK-N-SH", "Loucy", "NB4",
    "MM.1S", "MCF7", "OCI-LY3", "DOHH2", "SU-DHL-6", "SK-N-MC",
    "LNCaP clone FGC", "KMS-11", "A549", "22Rv1", "VCaP", "Jurkat",
]

headers = {"accept": "application/json"}

base_url = "https://www.ebi.ac.uk/biostudies/api/v1/search"


POLYA_KEYWORDS = ["poly-a", "polya", "oligo-dt"]
TOTAL_KEYWORDS = ["total rna", "ribo-depletion", "ribozero"]

TREATMENT_PROTOCOL_NAME_KEYWORDS = ["treatment protocol", "treatment "]
TREATMENT_TEXT_KEYWORDS = [
    "treated with", "treatment with", "vehicle control", "dmso",
    "knockdown", "knockout", "sirna", "shrna", "crispr",
    "overexpression", "transfected with", "transduced with",
]



def get_study_type(record_json):

    section = record_json.get("section", {})

    def scan_attributes(node):
        found = []
        for attr in node.get("attributes", []):
            if attr.get("name", "").strip().lower() == "study type":
                found.append(attr.get("value", ""))
        for subsection in node.get("subsections", []):
            if isinstance(subsection, dict):
                found.extend(scan_attributes(subsection))
        return found

    values = scan_attributes(section)
    return values  # list of study type strings found


def get_protocol_names(record_json):

    section = record_json.get("section", {})

    def scan(node):
        found = []
        node_type = node.get("type", "")
        if isinstance(node_type, str):
            found.append(node_type)
        for attr in node.get("attributes", []):
            if attr.get("name", "").strip().lower() == "type":
                found.append(attr.get("value", ""))
        for subsection in node.get("subsections", []):
            if isinstance(subsection, dict):
                found.extend(scan(subsection))
        return found

    return [v.strip().lower() for v in scan(section) if isinstance(v, str)]


def has_treatment_signal(record_json, full_text_lower):

    protocol_names = get_protocol_names(record_json)

    for name in protocol_names:
        if any(kw in name for kw in TREATMENT_PROTOCOL_NAME_KEYWORDS):
            return True, f"treatment protocol present: '{name}'"

    for kw in TREATMENT_TEXT_KEYWORDS:
        if kw in full_text_lower:
            return True, f"treatment keyword found in record text: '{kw}'"

    return False, ""


def classify_rna_type(accession):

    record_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}"

    try:
        res = requests.get(record_url, headers=headers)
    except requests.RequestException as e:
        return None, f"request error: {e}"

    if res.status_code != 200:
        return None, f"HTTP {res.status_code}"

    try:
        record_json = res.json()
    except ValueError:
        return None, "could not parse record JSON"

    study_types = get_study_type(record_json)

    if not study_types:
        return None, "no 'Study type' attribute found in record"

    study_types_lower = [st.strip().lower() for st in study_types]

    is_sequencing = any("rna-seq" in st for st in study_types_lower)

    if not is_sequencing:
        return None, f"not RNA-seq by Study type: {study_types}"

    # confirmed sequencing-based RNA-seq study
    text = res.text.lower()

    is_treated, treatment_reason = has_treatment_signal(record_json, text)
    if is_treated:
        return None, f"excluded as treated: {treatment_reason}"

    has_polya = any(kw in text for kw in POLYA_KEYWORDS)
    has_total = any(kw in text for kw in TOTAL_KEYWORDS)

    if has_polya and has_total:
        return "ambiguous", "both polyA and total RNA keywords present"
    elif has_polya:
        return "polyA", ""
    elif has_total:
        return "total RNA", ""
    else:
        return None

all_hits = []
zero_hit_lines = []

for cl in missing_cell_lines:

    query_string = f'"{cl}" RNA-seq'  # quoted name for precision, as established previously

    res = requests.get(
        base_url,
        headers=headers,
        params={"query": query_string, "type": "study"},
    )

    if res.status_code != 200:
        print(f"[{cl}] request failed: HTTP {res.status_code} - {res.text[:200]}")
        continue

    try:
        data = res.json()
    except ValueError:
        print(f"[{cl}] could not parse response body as JSON")
        continue

    hits = data.get("hits", [])

    if len(hits) == 0:
        zero_hit_lines.append(cl)
        print(f"[{cl}] 0 hits")
        time.sleep(0.5)
        continue

    print(f"[{cl}] {len(hits)} hits")

    for h in hits:

        acc = h.get("accession", "")

        rna_type, rna_type_note = classify_rna_type(acc)
        time.sleep(0.5)  # be polite - one extra request per hit now

        if rna_type is None:
            print(f"  [{acc}] skipping ({rna_type_note})")
            continue

        all_hits.append({
            "queried_cell_line": cl,
            "accession": acc,
            "title": h.get("title", ""),
            "release_date": h.get("release_date", ""),
            "rna_type": rna_type,
            "rna_type_note": rna_type_note,
        })

    time.sleep(0.5)  # be polite to the API

df = pd.DataFrame(all_hits)



df.to_csv("arrayexpress_missing_cell_lines_rnaseq.tsv", sep="\t", index=False)

print(df)