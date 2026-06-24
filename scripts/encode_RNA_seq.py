

import time, requests, pandas as pd

cell_lines = [
    "SJSA1","MG63","HL-60","NCI-H929","HepG2","HCT116","Panc1","DND-41",
    "SK-N-SH","K562","Loucy","NB4","PC-3","Karpas-422","MM.1S","MCF7",
    "OCI-LY3","DOHH2","A673","SU-DHL-6","SK-N-MC","LNCaP clone FGC",
    "KMS-11","A549","22Rv1","VCaP","Jurkat"
]

hits = []

for cl in cell_lines:

    r = requests.get(
        "https://www.encodeproject.org/search/",
        params={
            "type":"Experiment",
            "assay_term_name":"RNA-seq",
            "biosample_ontology.term_name":cl,
            "status":"released",
            "format":"json",
            "limit":"all"
        },
        headers={"accept":"application/json"}
    )

    for exp in r.json().get("@graph", []):

        acc = exp["accession"]

        d = requests.get(
            f"https://www.encodeproject.org/experiments/{acc}/?format=json",
            headers={"accept":"application/json"}
        ).json()

        reps = d.get("replicates", [])
        if not reps:
            continue

        untreated = all(
            not bs.get("treatments") and not bs.get("applied_modifications")
            for bs in (
                rep.get("library", {}).get("biosample", {})
                for rep in reps
            )
        )

        whole_cell = all(
            bs.get("subcellular_fraction_term_name") in [None, "whole cell", "Whole cell"]
            for bs in (
                rep.get("library", {}).get("biosample", {})
                for rep in reps
            )
        )

        if untreated and whole_cell:
            hits.append({
                "cell_line": cl,
                "accession": acc,
                "rna_type": d.get("assay_title"),
                "lab": d.get("lab", {}).get("title")
            })

        time.sleep(0.2)

df = pd.DataFrame(hits)
df.to_csv("encode_untreated_wholecell_rnaseq.tsv", sep="\t", index=False)
print(df)