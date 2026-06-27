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
    print(cl, "search status:", r.status_code)
    if r.status_code == 404:
        # valid empty result for this cell line, not an error
        continue
    exps = r.json().get("@graph", [])

    for exp in exps:
        acc = exp["accession"]
        d = requests.get(
            f"https://www.encodeproject.org/experiments/{acc}/?format=json",
            headers={"accept":"application/json"}
        ).json()
        time.sleep(0.2)

        reps = d.get("replicates", [])
        if not reps:
            continue

        biosamples = [
            rep.get("library", {}).get("biosample", {})
            for rep in reps
        ]

        untreated = all(
            not bs.get("treatments") and not bs.get("applied_modifications")
            for bs in biosamples
        )
        whole_cell = all(
            bs.get("subcellular_fraction_term_name") in [None, "whole cell", "Whole cell"]
            for bs in biosamples
        )

        if not (untreated and whole_cell):
            continue

        lab_name = d.get("lab", {}).get("title")

        # lab filter: only Thomas Gingeras (CSHL) and Barbara Wold (Caltech)
        if lab_name not in ["Thomas Gingeras, CSHL", "Barbara Wold, Caltech"]:
            continue

        # pull all files attached to this experiment
        files = d.get("files", [])

        if not files:
            print("no files found for", acc)
            continue

        for f in files:
            output_type = f.get("output_type")
            if output_type != "gene quantifications":
                continue  # skip BAMs, bigWigs, transcript-level files, etc

            file_acc = f.get("accession")
            assembly = f.get("assembly")
            genome_annotation = f.get("genome_annotation")
            file_format = f.get("file_format")
            preferred_default = f.get("preferred_default", False)
            href = f.get("href")
            file_url = f"https://www.encodeproject.org{href}" if href else None

            # biological replicate number(s) this file corresponds to
            bio_reps = f.get("biological_replicates", [])
            bio_rep_str = ",".join(str(x) for x in bio_reps)

            hits.append({
                "cell_line": cl,
                "accession": acc,
                "rna_type": d.get("assay_title"),
                "lab": lab_name,
                "file_accession": file_acc,
                "file_format": file_format,
                "output_type": output_type,
                "assembly": assembly,
                "genome_annotation": genome_annotation,
                "biological_replicates": bio_rep_str,
                "preferred_default": preferred_default,
                "file_url": file_url
            })

    time.sleep(0.2)

df = pd.DataFrame(hits)
df.to_csv("encode_untreated_wholecell_rnaseq_with_files.tsv", sep="\t", index=False)
print(df)
