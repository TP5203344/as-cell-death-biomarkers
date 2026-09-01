# -*- coding: utf-8 -*-
"""
Drug prediction for immune subtypes (K3) using Enrichr API
- Input: 17 core cell death genes (signature from LASSO + RF-RFE)
- Libraries: DSigDB, HDSigDB_Human_2021, LINCS_L1000_Chem_Pert_up/down
- Method: Enrichr enrichment (Fisher exact + z-score + combined score)

Requirements: requests, pandas
"""

import time, csv
import requests
import pandas as pd

ENRICHR = "https://maayanlab.cloud/Enrichr"
HDRS = {"User-Agent": "Mozilla/5.0"}


def add_list(genes):
    """Submit gene list (must use multipart/form-data, urlencoded returns 400)"""
    r = requests.post(f"{ENRICHR}/addList",
                      files={"list": (None, "\n".join(genes))},
                      headers=HDRS, timeout=120)
    r.raise_for_status()
    return r.json()["userListId"]


def enrich(uid, lib, max_tries=8):
    """Query enrichment results, waiting for server processing"""
    for attempt in range(max_tries):
        try:
            r = requests.get(f"{ENRICHR}/enrich",
                             params={"userListId": uid, "backgroundType": lib},
                             headers=HDRS, timeout=180)
            r.raise_for_status()
            data = r.json()
            total = sum(len(v) for v in data.values())
            if total > 0:
                return data
            print(f"    [{lib}] waiting ({attempt + 1}/{max_tries})...")
        except Exception as e:
            print(f"    [{lib}] retry ({attempt + 1}): {e}")
        time.sleep(10)
    return {}


def run_drug_prediction(core_genes, libs, out_csv):
    all_rows = []
    uid = add_list(core_genes)
    print(f"userListId={uid}")
    time.sleep(3)
    for lib in libs:
        res = enrich(uid, lib)
        for libname, arr in res.items():
            for row in arr:
                if len(row) >= 7:
                    all_rows.append({
                        "library": libname,
                        "term": row[1], "p": row[2], "z": row[3],
                        "combined": row[4], "genes": row[5], "adj_p": row[6],
                        "n_input": len(core_genes),
                    })
        n = sum(1 for r in all_rows if r["library"] == lib)
        print(f"  [{lib}] {n} results")
        time.sleep(2)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["library", "term", "p", "adj_p", "z_score", "combined_score", "genes", "n_input"])
        for r in sorted(all_rows, key=lambda x: x["p"]):
            w.writerow([r["library"], r["term"], f"{r['p']:.6e}",
                        f"{r['adj_p']:.6e}", f"{r['z']:.4f}", f"{r['combined']:.4f}",
                        r["genes"], r["n_input"]])
    return all_rows


if __name__ == "__main__":
    # 17 core cell death genes (from LASSO + RF-RFE, Figure 3)
    core17 = ["AIM2", "BMF", "CASP1", "CTSD", "CTSG", "CYBA", "IL18", "ITGA5",
              "ITGAM", "MPO", "NCF2", "NLRP3", "PGD", "PRKAA2", "PYCARD",
              "RIPK3", "ZEB1"]
    libs = ["DSigDB", "HDSigDB_Human_2021",
            "LINCS_L1000_Chem_Pert_up", "LINCS_L1000_Chem_Pert_down"]
    rows = run_drug_prediction(core17, libs, "drug_enrichr_results.csv")
    sig = [r for r in rows if r["adj_p"] < 0.05]
    print(f"Significant (adjP<0.05): {len(sig)}")
