import re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt


COMMERCIAL = ["Erica", "Kadie", "Malgas", "Nemesia"]
SA = COMMERCIAL + ["UVC8", "RGT_Planet"] 
INTL = ["Barke", "Golden_Promise", "Hocket", "Igri", "Maximus", "Morex"]
ORDER = ["Morex", "Maximus", "Igri", "Barke", "Golden_Promise", "Hocket",
         "RGT_Planet", "UVC8", "Erica", "Malgas", "Kadie", "Nemesia"]
COL = {"Core": "#1f77b4", "Shell": "#ff7f0e", "Cloud": "#2ca02c"}

# loads file
df = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
if df.shape[0] > df.shape[1]: df = df.T
df.index = [re.sub(r"\.fna(\.nucl)?$", "", i) for i in df.index]
df.columns = [re.sub(r"\.fna$", "", re.sub(r"^\d+_", "", c)) for c in df.columns]
P = df.apply(pd.to_numeric, errors="coerce").fillna(0).gt(0)
N, total = P.shape

# core, shell, cloud
n = P.sum()
cat = pd.Series("Shell", index=P.columns)   # shell = 2 to N-1 accessions
cat[n == 1] = "Cloud"                       # cloud = 1 accession
cat[n == N] = "Core"                        # core  = all N accessions
log = [f"Total clusters: {total:,}"]
log += [f"{c}: {(cat == c).sum():,} ({100 * (cat == c).mean():.1f}%)" for c in COL]

per_acc = pd.DataFrame({c: P.loc[:, cat == c].sum(axis=1) for c in COL})
log.append("\nClusters per accession:\n" + per_acc.assign(Total=per_acc.sum(axis=1)).to_string())

# chromosome specific grouping
cmap = dict(pd.read_csv(sys.argv[2], sep="\t", header=None).values) if len(sys.argv) > 2 else {}
def chrom(name):
    if name in cmap: return cmap[name]
    m = re.search(r"\.([1-7])HG|chr([1-7])H", name)
    return (m.group(1) or m.group(2)) + "H" if m else "Unknown"
chrs = pd.Series([chrom(c) for c in P.columns], index=P.columns)
CHR = [f"{i}H" for i in range(1, 8)] + ["Unknown"]
per_chr = pd.DataFrame({c: cat[cat == c].groupby(chrs).size() for c in COL}).reindex(CHR).fillna(0).astype(int)

# Group-specific clusters
sa_only = P[SA].any() & ~P[INTL].any()
c5 = (chrs[sa_only] == "5H").sum()
log.append(f"\nSA-specific clusters (>=1 SA, no international): {sa_only.sum():,}; "
           f"on 5H: {c5} ({100 * c5 / sa_only.sum():.1f}%)")
rgt = P["RGT_Planet"] & ~P[INTL].any()
log.append(f"RGT_Planet-specific vs international: {rgt.sum():,}")
no_uvc8 = P[COMMERCIAL].any() & ~P["UVC8"]
all4 = no_uvc8 & P[COMMERCIAL].all()
log.append(f"Commercial SA clusters absent in UVC8: {no_uvc8.sum():,}; in all four: "
           f"{all4.sum():,} ({100 * all4.sum() / no_uvc8.sum():.1f}%)")

# Jaccard
X = P.values.astype(int)
inter = X @ X.T
s = X.sum(1)
J = pd.DataFrame(inter / (s[:, None] + s[None, :] - inter), index=P.index, columns=P.index)
def rng(a, b):
    v = [J.loc[x, y] for x in a for y in b if x != y]
    return f"{min(v):.2f}-{max(v):.2f}"
log.append("\nJaccard ranges:")
log.append(f"  SA (commercial + UVC8) vs international: {rng(COMMERCIAL + ['UVC8'], INTL)}")
log.append(f"  commercial vs commercial: {rng(COMMERCIAL, COMMERCIAL)}")
log.append(f"  UVC8 vs commercial: {rng(['UVC8'], COMMERCIAL)}")
log.append(f"  RGT_Planet vs international: {rng(['RGT_Planet'], INTL)}")
log.append(f"  RGT_Planet vs commercial: {rng(['RGT_Planet'], COMMERCIAL)}")

# Accumulation curve
rs = np.random.default_rng(1)
orders = [rs.permutation(N) for _ in range(200)]
curve = np.mean([[P.values[o[:k]].any(axis=0).sum() for k in range(1, N + 1)] for o in orders], axis=0)
log.append("\nAccumulation: " + ", ".join(f"{k}:{v:,.0f}" for k, v in enumerate(curve, 1)))
log.append(f"Clusters added by last accession: {curve[-1] - curve[-2]:,.0f}")

open("summary_stats.txt", "w").write("\n".join(log))
print("\n".join(log))

# Figures
def stacked(ax, d, title, xlabel):
    bot = np.zeros(len(d))
    for c in COL:
        ax.bar(d.index, d[c], bottom=bot, color=COL[c], width=0.55, label=c)
        for i, (v, b) in enumerate(zip(d[c], bot)):
            if v > 0.03 * d.sum(axis=1).max():
                ax.text(i, b + v / 2, f"{v:,}", ha="center", va="center", fontsize=6, color="white")
        bot += d[c].values
    ax.set(title=title, xlabel=xlabel, ylabel="Number of ortholog groups")
    ax.legend(fontsize=7, loc="lower left")

fig, ax = plt.subplots(figsize=(8, 5))
stacked(ax, per_acc.sort_index(), "Core, Shell and Cloud Gene Content", "Accession")
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
fig.tight_layout(); fig.savefig("fig3A_per_accession.png", dpi=300)

fig, ax = plt.subplots(figsize=(7, 5))
stacked(ax, per_chr, "Core, Shell and Cloud Gene Content per Chromosome", "Chromosome")
fig.tight_layout(); fig.savefig("fig3B_per_chromosome.png", dpi=300)

M = J.loc[ORDER, ORDER]
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(M, cmap="viridis", vmax=1)
ax.set_xticks(range(N)); ax.set_yticks(range(N))
ax.set_xticklabels(ORDER, rotation=45, ha="right"); ax.set_yticklabels(ORDER)
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_color("teal" if lab.get_text() in SA else "orangered")
for i in range(N):
    for j in range(N):
        ax.text(j, i, f"{M.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7,
                color="black" if M.iloc[i, j] > 0.9 else "white")
ax.set_title("Gene-content similarity\n(teal = SA-grown, orange = international)")
fig.colorbar(im, fraction=0.046).set_label("Jaccard similarity")
fig.tight_layout(); fig.savefig("fig4A_jaccard.png", dpi=300)

fig, ax = plt.subplots(figsize=(5.5, 4.5))
ax.plot(range(1, N + 1), curve, "o-")
ax.set(xticks=range(1, N + 1), xlabel="Number of accessions", ylabel="Number of clusters (pan-genome)")
fig.tight_layout(); fig.savefig("fig4B_accumulation.png", dpi=300)
