# -*- coding: utf-8 -*-
"""GSE163154 (MaasHPS) 17 核心基因外部验证 - 修正版"""
import sys, os
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

datadir = r"C:\Users\19267\.qclaw\workspace-agent-b69b3ebd\as_bioinfo\GSE163154_data"
sm = os.path.join(datadir, "GSE163154_series_matrix.txt")
mapf = os.path.join(datadir, "GPL6104_ilmn2symbol.txt")

# 1. 读取映射
ilmn2sym = {}
with open(mapf, encoding="utf-8") as f:
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) == 2 and p[0] and p[1]:
            ilmn2sym[p[0]] = p[1].strip().upper()

core17 = ["AIM2","BMF","CASP1","CTSD","CTSG","CYBA","IL18","ITGA5","ITGAM",
          "MPO","NCF2","NLRP3","PGD","PRKAA2","PYCARD","RIPK3","ZEB1"]

with open(sm, encoding="utf-8", errors="ignore") as f:
    lines = f.readlines()

# 2. 样本分组：从所有 characteristics 行中找 classification
group = None
for ln in lines:
    if ln.startswith('!Sample_characteristics_ch1'):
        vals = [t.strip('"') for t in ln.strip().split("\t")[1:]]
        if vals and any("classification" in v or "IPH" in v or "non-IPH" in v for v in vals):
            group = [1 if "IPH" in v and "non" not in v else 0 for v in vals]
            print(f"分组行: {vals[:5]}...")

if group is None:
    print("未找到分组行!"); sys.exit(1)

print(f"样本数: {len(group)}, IPH: {sum(group)}, non-IPH: {len(group)-sum(group)}")

# 3. 解析矩阵（跳过表头行）
in_table = False
expr_rows = {}
n = 0
for ln in lines:
    if 'series_matrix_table_begin' in ln:
        in_table = True
        continue
    if 'series_matrix_table_end' in ln:
        break
    if in_table:
        parts = ln.rstrip("\n").split("\t")
        if parts[0].startswith('"ID_REF"') or parts[0] == "ID_REF":
            continue  # 表头
        probe = parts[0].strip('"')
        try:
            vals = [float(v.strip('"')) for v in parts[1:]]
        except ValueError:
            continue
        sym = ilmn2sym.get(probe, "")
        if sym in core17:
            expr_rows.setdefault(sym, []).append(vals)
        n += 1

print(f"总行数: {n}, 命中基因数: {len(expr_rows)}")

# 多探针取均值
expr = {}
for sym, vlist in expr_rows.items():
    arr = np.array(vlist)
    expr[sym] = arr.mean(axis=0)

# 4. 差异检验
import scipy.stats as st
print(f"\n{'Gene':<10}{'IPH_mean':>10}{'nonIPH_mean':>12}{'delta':>10}{'P':>12}")
results = []
for sym in core17:
    if sym not in expr:
        print(f"{sym:<10} 缺失")
        continue
    v = expr[sym]
    iph = v[np.array(group) == 1]
    non = v[np.array(group) == 0]
    d = iph.mean() - non.mean()
    t, p = st.ttest_ind(iph, non)
    results.append((sym, d, p))
    print(f"{sym:<10}{iph.mean():>10.4f}{non.mean():>12.4f}{d:>10.4f}{p:>12.4e}")

# 5. ROC
from sklearn.metrics import roc_auc_score
if len(results) >= 5:
    score = np.mean([expr[s] for s, _, _ in results], axis=0)
    auc = roc_auc_score(group, score)
    print(f"\n17基因均值评分 AUC (IPH识别): {auc:.4f}")
    import json
    out = {"n_samples": len(group), "n_IPH": sum(group), "n_nonIPH": len(group)-sum(group),
           "auc_17gene_mean": auc,
           "per_gene": [{"gene": s, "delta": round(d, 4), "p": p} for s, d, p in results]}
    with open(os.path.join(datadir, "GSE163154_validation_results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("已保存 GSE163154_validation_results.json")
