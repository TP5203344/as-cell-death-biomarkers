# -*- coding: utf-8 -*-
"""生成 Figure 6 D 面板: 异质性(Cochran's Q) + 多效性(MR-Egger intercept) 可视化"""
import csv, os, sys
sys.stdout.reconfigure(encoding="utf-8")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

mr = r"C:\Users\19267\Desktop\纯生信文档\0.顶刊复刻代码\Figure 4 免疫浸润与诊断模型\MR"
outdir = r"C:\Users\19267\.qclaw\workspace-agent-b69b3ebd\as_bioinfo"

genes = ["AIM2","BMF","CASP1","CTSD","CTSG","CYBA","ITGA5","ITGAM","MPO","NLRP3"]

# 读取数据
data = []
for g in genes:
    hf = os.path.join(mr, f"mr_heterogeneity_{g}.csv")
    pf = os.path.join(mr, f"mr_pleiotropy_{g}.csv")
    hrows = list(csv.DictReader(open(hf, encoding="utf-8")))
    prows = list(csv.DictReader(open(pf, encoding="utf-8")))
    ivw = [r for r in hrows if "Inverse variance" in r.get("method","")]
    pleio = prows[0] if prows else {}
    q_pval = float(ivw[0]["Q_pval"]) if ivw else np.nan
    intercept = float(pleio.get("egger_intercept", np.nan))
    se = float(pleio.get("se", np.nan))
    pval = float(pleio.get("pval", np.nan))
    data.append(dict(gene=g, q_pval=q_pval, intercept=intercept, se=se, pval=pval))

# 按基因排序（保持输入顺序，可读性）
genes_plot = [d["gene"] for d in data]
y = np.arange(len(genes_plot))[::-1]  # 从上到下

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300,
                                gridspec_kw={"width_ratios":[1.15, 1]})

# ========== 左: 多效性森林图 (Egger intercept ± 95% CI) ==========
for i, d in enumerate(data):
    color = "#d62728" if d["pval"] < 0.05 else "#1f77b4"
    ax1.errorbar(d["intercept"], y[i], xerr=1.96*d["se"], fmt="o",
                 color=color, ecolor=color, capsize=3, markersize=5, elinewidth=1.2)
    # P 值标注
    pstr = f"P={d['pval']:.3f}"
    ax1.text(max(d["intercept"]+1.96*d["se"], 0.001) + 0.008, y[i]+0.18,
             pstr, fontsize=7, va="center", color="0.3")

ax1.axvline(0, color="0.4", linestyle="--", linewidth=1)
ax1.set_yticks(y)
ax1.set_yticklabels(genes_plot, fontsize=9)
ax1.set_xlabel("MR-Egger intercept (pleiotropy)", fontsize=10)
ax1.set_title("Pleiotropy test (MR-Egger intercept)", fontsize=11, fontweight="bold")
ax1.tick_params(axis="x", labelsize=8)
ax1.grid(axis="x", alpha=0.25, linestyle=":")

# ========== 右: 异质性条形图 (-log10 IVW Q P) ==========
neglog = [-np.log10(d["q_pval"]) if d["q_pval"] > 0 else 30 for d in data]
colors = ["#d62728" if d["q_pval"] < 0.05 else "#1f77b4" for d in data]
ax2.barh(y, neglog, color=colors, alpha=0.85, edgecolor="none")
ax2.axvline(-np.log10(0.05), color="0.4", linestyle="--", linewidth=1)
ax2.text(-np.log10(0.05)+0.15, -0.5, "P = 0.05", fontsize=8, color="0.35")
ax2.set_yticks(y)
ax2.set_yticklabels(genes_plot, fontsize=9)
ax2.set_xlabel("-log10(P) of Cochran's Q (IVW)", fontsize=10)
ax2.set_title("Heterogeneity test (Cochran's Q)", fontsize=11, fontweight="bold")
ax2.tick_params(axis="x", labelsize=8)
ax2.grid(axis="x", alpha=0.25, linestyle=":")
# 标注显著值
for i, d in enumerate(data):
    if d["q_pval"] < 0.05:
        ax2.text(neglog[i]+0.15, y[i], f"P={d['q_pval']:.1e}", fontsize=7, va="center")

plt.tight_layout()

svg_path = os.path.join(outdir, "MR_het_pleio_panel.svg")
png_path = os.path.join(outdir, "MR_het_pleio_panel.png")
plt.savefig(svg_path, format="svg", bbox_inches="tight")
plt.savefig(png_path, format="png", bbox_inches="tight", dpi=300)
plt.close()

print("D 面板已生成:")
print("  SVG:", svg_path, os.path.getsize(svg_path), "B")
print("  PNG:", png_path, os.path.getsize(png_path), "B")

# 打印数据确认
print("\n=== 数据确认 ===")
for d in data:
    sig_het = "HET*" if d["q_pval"] < 0.05 else ""
    sig_pleio = "PLEIO*" if d["pval"] < 0.05 else ""
    print(f"{d['gene']:6s} Q_P={d['q_pval']:.2e} {sig_het:6s} intercept={d['intercept']:+.3f} P={d['pval']:.3f} {sig_pleio}")
