# -*- coding: utf-8 -*-
"""
方案A Step③: 免疫浸润分析 (CIBERSORT nnls 实现)
- 使用 GSE100927 全量表达矩阵 (34964 基因, 69 AS vs 35 Control)
- LM22 免疫细胞签名反卷积
- 输出: 免疫细胞比例、AS vs Control 差异、核心基因-免疫细胞相关性
"""
import pandas as pd, numpy as np, sys, os, json, subprocess
from scipy.optimize import nnls
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding='utf-8')

BASE = r"C:\Users\19267\Desktop\纯生信文档\0.顶刊复刻代码"
OUTDIR = os.path.join(BASE, "Figure 4 免疫浸润与诊断模型")
os.makedirs(OUTDIR, exist_ok=True)

# ============ LM22 签名 (基因 -> 细胞类型) ============
lm22_signatures = {
    "B cells naive": ["CD79A","MS4A1","CD19","CD79B","BLK","FCRL1","TCL1A","VPREB3","SPIB","BANK1"],
    "B cells memory": ["CD27","CD40","CD80","CR2","TNFRSF13B","TNFRSF13C","TNFRSF17","LTB","POU2AF1","IRF4"],
    "Plasma cells": ["XBP1","SDC1","TNFRSF17","MZB1","JCHAIN","FKBP11","DERL3","SEC61B","TNFRSF13B","SSR4"],
    "T cells CD8": ["CD8A","CD8B","GZMA","GZMB","GZMK","PRF1","IFNG","CCL5","CXCR3","NKG7"],
    "T cells CD4 naive": ["CCR7","LEF1","TCF7","SELL","IL7R","MAL","CD27","FOXP1","TRAT1","RGCC"],
    "T cells CD4 memory resting": ["CD28","IL7R","CCR7","CD27","SELL","TCF7","LEF1","MAL","S100A4","KLRB1"],
    "T cells CD4 memory activated": ["CD69","ICOS","CD40LG","IFNG","TNF","IL2","CXCR3","CCR5","IL4","IL13"],
    "T cells follicular helper": ["CXCR5","PDCD1","BCL6","ICOS","IL21","MAF","TNFRSF4","CD40LG","CXCL13","SH2D1A"],
    "T cells regulatory (Tregs)": ["FOXP3","IL2RA","CTLA4","IKZF2","IKZF4","TNFRSF18","CCR8","BATF","ENTPD1","LAX1"],
    "T cells gamma delta": ["TRDC","TRGC1","TRGC2","TRDV1","KLRD1","NKG7","KLRC1","CD160","FCRL6","CCL4"],
    "NK cells resting": ["KLRF1","KLRD1","NCR1","KLRB1","KLRC1","KLRG1","SH2D1B","FCGR3A","NKG7","GZMB"],
    "NK cells activated": ["IFNG","GZMB","PRF1","KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2","FCGR3A","NKG7","CX3CR1"],
    "Monocytes": ["CD14","LYZ","S100A8","S100A9","FCGR3A","CSF1R","CD68","ITGAM","CCR2","TLR4"],
    "Macrophages M0": ["CD68","CD14","FCGR1A","FCGR1B","CSF1R","ITGAM","LYZ","MARCO","MRC1","MSR1"],
    "Macrophages M1": ["NOS2","IL1B","TNF","IL6","CXCL9","CXCL10","CXCL11","CCL5","CD80","CD86","HLA-DRA","HLA-DRB1","IL12A","IL12B","IRF5","IRF1","TLR2","TLR4","STAT1","GBP2"],
    "Macrophages M2": ["MRC1","CD163","MSR1","IL10","TGFB1","ARG1","CCL18","CCL22","CD209","CXCL13","FN1","SPP1","MMP9","MMP12","PPARG","IRF4","STAT6","CCL17","CCL13","CCL23"],
    "Dendritic cells resting": ["CD1C","CD1A","CD1B","CD1D","CLEC10A","FCER1A","FLT3","ITGAX","CD207","CD209"],
    "Dendritic cells activated": ["CCL19","CCL21","CCR7","CD40","CD80","CD83","CD86","LAMP3","BATF3","IRF8"],
    "Mast cells resting": ["TPSAB1","TPSB2","CPA3","MS4A2","HDC","TPSD1","CMA1","CTSG","GATA2","KIT"],
    "Mast cells activated": ["IL13","IL4","IL5","IL6","IL8","CCL2","CCL3","CCL4","TNF","VEGFA"],
    "Eosinophils": ["CLC","PRG2","PRG3","EPX","RNASE3","CCR3","IL5RA","SIGLEC8","OLFM4","ALOX15"],
    "Neutrophils": ["FCGR3B","CSF3R","S100A8","S100A9","S100A12","CEACAM8","ITGAM","ITGB2","CXCR1","CXCR2","FPR1","FPR2","MMP8","MMP9","MPO","ELANE","LTF","LCN2","ARG1","CD177"],
}

# 反转: gene -> [cell types]
gene2cells = {}
for cell, genes in lm22_signatures.items():
    for g in genes:
        gene2cells.setdefault(g, []).append(cell)
all_sig_genes = sorted(gene2cells.keys())

# ============ 加载数据 ============
expr = pd.read_csv(r"C:\Users\19267\.qclaw\workspace-agent-b69b3ebd\as_bioinfo\_gse100927_expr_full.csv", index_col=0)
print(f"GSE100927 表达矩阵: {expr.shape}")

# 探针 -> Symbol 映射 (gene_map)
import json
r_script_map = '''
x <- readRDS("C:/Users/19267/Desktop/纯生信文档/0.顶刊复刻代码/Figure 1 研究设计与转录组差异分析/热图GSE100927/expr_norm.rds")
gm <- x[["gene_map"]]
df <- data.frame(probe=names(gm), symbol=unname(gm), stringsAsFactors=FALSE)
write.csv(df, "C:/Users/19267/.qclaw/workspace-agent-b69b3ebd/as_bioinfo/_probe_symbol_map.csv", row.names=FALSE, quote=FALSE)
'''
proc = subprocess.run(["Rscript", "-e", r_script_map], capture_output=True, text=True, encoding='utf-8', errors='replace')
if proc.returncode != 0:
    print("STDERR:", proc.stderr[:1500])
    sys.exit(1)
ps_map = pd.read_csv(r"C:\Users\19267\.qclaw\workspace-agent-b69b3ebd\as_bioinfo\_probe_symbol_map.csv")
probe2sym = dict(zip(ps_map['probe'], ps_map['symbol']))

# 转换到 Symbol (保留表达矩阵中存在的探针)
expr_sym = expr.copy()
expr_sym.index = [probe2sym.get(p, p) for p in expr.index]
expr_sym = expr_sym[~expr_sym.index.duplicated(keep='first')]
# 去除非标准符号 (含点/ENST等)
mask = [isinstance(s, str) and s == s.upper() and not s.startswith('ENST') and not s.startswith('LOC') and not s.startswith('XLOC') for s in expr_sym.index]
expr_sym = expr_sym[mask]
print(f"转换后表达矩阵: {expr_sym.shape}")

group = pd.read_csv(r"C:\Users\19267\.qclaw\workspace-agent-b69b3ebd\as_bioinfo\_gse100927_group.csv")
group.columns = ['sample', 'title', 'group']
group['group01'] = (group['group'] == 'Atherosclerotic').astype(int)
group_map = dict(zip(group['sample'], group['group']))

# 签名基因命中检查
common = [g for g in all_sig_genes if g in expr_sym.index]
print(f"签名基因 {len(all_sig_genes)} 个, 表达矩阵命中 {len(common)} 个")

# ============ CIBERSORT 反卷积 ============
def cibersort(expr_df, gene2cells_map, common_genes, n_perm=100, seed=123):
    rng = np.random.RandomState(seed)
    cells = [c for c in lm22_signatures.keys() if any(g in common_genes for g in lm22_signatures[c])]
    sig_mat = np.zeros((len(common_genes), len(cells)))
    for j, cell in enumerate(cells):
        for i, g in enumerate(common_genes):
            if cell in gene2cells_map[g]:
                sig_mat[i, j] = 1.0
    # 每个细胞至少 1 个基因
    keep = sig_mat.sum(axis=0) > 0
    cells = [c for c, k in zip(cells, keep) if k]
    sig_mat = sig_mat[:, keep]
    
    X = expr_df.loc[common_genes].values.astype(float)
    n_samples = X.shape[1]
    fracs = np.zeros((n_samples, len(cells)))
    pvals = np.zeros(n_samples)
    
    for s in range(n_samples):
        x = X[:, s].copy()
        x = x - x.min() + 1e-6
        x = x / x.sum()
        coef, _ = nnls(sig_mat, x)
        total = coef.sum()
        if total > 0:
            coef = coef / total
        fracs[s, :] = coef
        if n_perm > 0:
            obs = np.corrcoef(sig_mat @ coef, x)[0, 1]
            cnt = 0
            for _ in range(n_perm):
                xp = rng.permutation(x)
                cp, _ = nnls(sig_mat, xp)
                csum = cp.sum()
                if csum > 0:
                    cp = cp / csum
                if np.corrcoef(sig_mat @ cp, xp)[0, 1] >= obs:
                    cnt += 1
            pvals[s] = (cnt + 1) / (n_perm + 1)
    
    out = pd.DataFrame(fracs, columns=cells, index=expr_df.columns)
    out['P-value'] = pvals
    return out

print("\n=== 免疫浸润反卷积 (GSE100927) ===")
imm = cibersort(expr_sym, gene2cells, common, n_perm=100)
print(f"完成: {imm.shape}")
imm.to_csv(os.path.join(OUTDIR, "GSE100927_immune_infiltration.csv"))
print(f"显著 (P<0.05): {(imm['P-value'] < 0.05).sum()} / {len(imm)} 样本")

# ============ 堆叠柱状图 ============
imm_plot = imm.drop(columns=['P-value']).T
fig, ax = plt.subplots(figsize=(16, 5))
bottom = np.zeros(imm_plot.shape[1])
colors = plt.cm.tab20(np.linspace(0, 1, imm_plot.shape[0]))
for i, (cell, row) in enumerate(imm_plot.iterrows()):
    ax.bar(range(imm_plot.shape[1]), row.values, bottom=bottom, label=cell, color=colors[i], width=0.8)
    bottom += row.values
ax.set_xlabel("Samples", fontsize=11)
ax.set_ylabel("Relative fraction", fontsize=11)
ax.set_title("Immune Cell Infiltration - GSE100927 (AS vs Control)", fontsize=13)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)
plt.subplots_adjust(bottom=0.18, left=0.08)
plt.savefig(os.path.join(OUTDIR, "GSE100927_immune_stack.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTDIR, "GSE100927_immune_stack.svg"), bbox_inches='tight')
plt.close()
print("堆叠柱状图已保存")

# ============ AS vs Control 差异 ============
imm_t = imm.drop(columns=['P-value']).copy()
imm_t['group'] = [group_map.get(s, 'Unknown') for s in imm_t.index]
results = []
for cell in imm_t.columns[:-1]:
    as_vals = imm_t.loc[imm_t['group'] == 'Atherosclerotic', cell]
    ctrl = imm_t.loc[imm_t['group'] == 'Control', cell]
    if len(as_vals) > 2 and len(ctrl) > 2:
        stat, p = stats.mannwhitneyu(as_vals, ctrl, alternative='two-sided')
        results.append({'Cell': cell, 'AS_mean': as_vals.mean(), 'Control_mean': ctrl.mean(),
                        'FC': as_vals.mean() / (ctrl.mean() + 1e-6), 'P': p,
                        'log2FC': np.log2(as_vals.mean() / (ctrl.mean() + 1e-6))})
res_df = pd.DataFrame(results).sort_values('P')
res_df.to_csv(os.path.join(OUTDIR, "immune_diff_AS_vs_Control.csv"), index=False)
print("\n=== AS vs Control 差异免疫细胞 (Top 15) ===")
print(res_df.head(15).to_string(index=False))

# 差异箱线图
sig_cells = res_df[res_df['P'] < 0.05].sort_values('P')['Cell'].tolist()
if len(sig_cells) > 8:
    sig_cells = sig_cells[:8]
if sig_cells:
    ncol = 4
    nrow = (len(sig_cells) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3.5*nrow))
    axes = np.array(axes).flatten()
    for ax, cell in zip(axes, sig_cells):
        as_vals = imm_t.loc[imm_t['group'] == 'Atherosclerotic', cell]
        ctrl = imm_t.loc[imm_t['group'] == 'Control', cell]
        ax.boxplot([ctrl, as_vals], tick_labels=['Control', 'AS'], widths=0.5)
        ax.set_title(cell, fontsize=9)
        ax.set_ylabel("Fraction")
    for ax in axes[len(sig_cells):]:
        ax.axis('off')
    plt.suptitle("Differential Immune Infiltration (AS vs Control)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "immune_diff_boxplot.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTDIR, "immune_diff_boxplot.svg"), bbox_inches='tight')
    plt.close()
    print(f"\n差异箱线图已保存: {len(sig_cells)} 个显著细胞")

# ============ 核心基因 vs 免疫细胞相关性 ============
core_genes = ["XAF1", "GZMB", "RNASE1", "TMEM106A", "DOK6", "CDH19", "CATSPERB", "ASPA", "TTYH2", "SCAMP5"]
core_available = [g for g in core_genes if g in expr_sym.index]
print(f"\n核心基因在表达矩阵中命中: {len(core_available)}/10 → {core_available}")

if core_available:
    corr_rows = []
    for gene in core_available:
        gexpr = expr_sym.loc[gene].values.astype(float)
        for cell in imm.columns[:-1]:
            cexpr = imm[cell].values.astype(float)
            r, p = stats.spearmanr(gexpr, cexpr)
            corr_rows.append({'Gene': gene, 'Cell': cell, 'rho': r, 'P': p})
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(os.path.join(OUTDIR, "core_gene_immune_correlation.csv"), index=False)
    
    # 相关性热图
    pivot = corr_df.pivot(index='Gene', columns='Cell', values='rho')
    pivot_p = corr_df.pivot(index='Gene', columns='Cell', values='P')
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(pivot.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=10)
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            val = pivot.values[i, j]
            pv = pivot_p.values[i, j]
            star = '***' if pv < 0.001 else ('**' if pv < 0.01 else ('*' if pv < 0.05 else ''))
            ax.text(j, i, f'{val:.2f}{star}', ha='center', va='center', fontsize=7)
    plt.colorbar(im, ax=ax, label='Spearman rho')
    ax.set_title("Core Genes vs Immune Cell Correlation (GSE100927)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "core_gene_immune_corr_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTDIR, "core_gene_immune_corr_heatmap.svg"), bbox_inches='tight')
    plt.close()
    print("核心基因-免疫细胞相关性热图已保存")

print("\n✅ Step③ 免疫浸润分析全部完成")
print(f"输出目录: {OUTDIR}")
