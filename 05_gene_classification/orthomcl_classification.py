import argparse
import os
import sys
from collections import Counter, defaultdict

# ---------------------------------------------------------------- CONFIG ----
DEFAULT_GROUPS = "C:/HonoursProj/ortholog_clustering/groups.txt"
DEFAULT_SINGLETONS = "C:/HonoursProj/ortholog_clustering/singletons.txt"

# Taxon names are compared case-insensitively. Add aliases if your file
# spells a taxon differently (e.g. "RGT_Planet", "RGTPlanet", "Planet").
SA_TAXA = {
    "erica": "Erica",
    "kadie": "Kadie",
    "malgas": "Malgas",
    "nemesia": "Nemesia",
    "rgt_planet": "RGT_Planet",
    "rgtplanet": "RGT_Planet",
    "planet": "RGT_Planet",
    "uvc8": "UVC8",
}
SA5 = {"Erica", "Kadie", "Malgas", "Nemesia", "RGT_Planet"}  # SA cultivars excl. UVC8
SA6 = SA5 | {"UVC8"}                                          # all SA-grown
RGT = "RGT_Planet"

# Substrings in a gene ID that mark it as an international cultivar.
# Used to detect mislabelled members such as "Erica|mRNA:Erica_BOWMAN.3HG...".
INTL_ID_MARKERS = ["MOREX", "IGRI", "HOCKETT", "GOLDEN_PROMISE", "BARKE", "MAXIMUS"]  # Bowman added because it appears in your Barley_OG_1000 line
# --------------------------------------------------------------------------


def normalise(taxon):
    return SA_TAXA.get(taxon.lower(), taxon)


def looks_international(gene_id):
    g = gene_id.upper()
    return any(m in g for m in INTL_ID_MARKERS)


def parse_member(member, fix_labels, mismatches, gid):
    taxon, _, gene = member.partition("|")
    taxon = normalise(taxon)
    if taxon in SA6 and looks_international(gene):
        mismatches.append((gid, member))
        if fix_labels:
            taxon = "INTERNATIONAL_MISLABELLED"
    return taxon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("groups", nargs="?", default=DEFAULT_GROUPS)
    ap.add_argument("singletons", nargs="?", default=DEFAULT_SINGLETONS)
    ap.add_argument("--outdir", default="orthomcl_results")
    ap.add_argument("--fix-labels", action="store_true",
                    help="Treat SA-labelled members with an international gene ID "
                         "(e.g. Erica|..._BOWMAN) as international instead of SA.")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    categories = {
        "sa_only_all6": [],        # every taxon in SA6
        "sa5_only_no_uvc8": [],    # every taxon in SA5, UVC8 absent
        "rgt_no_international": [],# RGT_Planet present, every taxon in SA6
        "rgt_planet_only": [],     # RGT_Planet is the only taxon
        "uvc8_only": [],
    }
    group_lines = {}               # gid -> original member string (only kept for hits)
    total_groups = 0
    taxon_presence = Counter()     # groups containing each taxon
    mismatches = []
    bad_lines = 0

    with open(args.groups) as f, \
         open(os.path.join(args.outdir, "group_taxa_matrix.tsv"), "w") as matrix:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if ":" not in line:
                bad_lines += 1
                continue
            header, members = line.split(":", 1)
            gid = header.split("(")[0].strip()
            members = members.split()
            total_groups += 1

            tx = {parse_member(m, args.fix_labels, mismatches, gid) for m in members}
            for t in tx:
                taxon_presence[t] += 1

            hit = False
            if tx <= SA6:
                categories["sa_only_all6"].append(gid); hit = True
                if tx == {RGT}:
                    categories["rgt_planet_only"].append(gid)
                if tx == {"UVC8"}:
                    categories["uvc8_only"].append(gid)
                if RGT in tx:
                    categories["rgt_no_international"].append(gid)
                if "UVC8" not in tx:
                    categories["sa5_only_no_uvc8"].append(gid)
            if hit:
                group_lines[gid] = " ".join(members)
            matrix.write(f"{gid}\t{len(members)}\t{','.join(sorted(tx))}\n")

    # ---- write category files
    for name, gids in categories.items():
        with open(os.path.join(args.outdir, f"{name}.txt"), "w") as o:
            for g in gids:
                o.write(f"{g}: {group_lines[g]}\n")

    # ---- singletons (one gene per line, or whitespace-separated)
    sing_taxa = Counter()
    n_sing = 0
    with open(args.singletons) as f:
        for line in f:
            for item in line.split():
                taxon = normalise(item.split("|")[0])
                gene = item.partition("|")[2]
                if taxon in SA6 and looks_international(gene):
                    taxon = "INTERNATIONAL_MISLABELLED" if args.fix_labels else taxon
                sing_taxa[taxon] += 1
                n_sing += 1

    # ---- report
    lines = []
    lines.append(f"Total groups:                                   {total_groups}")
    lines.append(f"Total singletons:                               {n_sing}")
    lines.append("")
    lines.append(f"SA-only groups (all 6 incl. UVC8):              {len(categories['sa_only_all6'])}")
    lines.append(f"Unique to 5 SA cultivars (no UVC8):             {len(categories['sa5_only_no_uvc8'])}")
    lines.append(f"Contain RGT_Planet, no international:           {len(categories['rgt_no_international'])}")
    lines.append(f"   of which RGT_Planet only:                    {len(categories['rgt_planet_only'])}")
    lines.append(f"UVC8-only groups:                               {len(categories['uvc8_only'])}")
    lines.append("")
    lines.append("Groups containing each taxon:")
    for t, c in sorted(taxon_presence.items(), key=lambda x: -x[1]):
        lines.append(f"   {t:28s}{c}")
    lines.append("")
    lines.append("Singletons per taxon:")
    for t, c in sorted(sing_taxa.items(), key=lambda x: -x[1]):
        lines.append(f"   {t:28s}{c}")
    sa_sing = sum(c for t, c in sing_taxa.items() if t in SA6)
    lines.append(f"   SA-grown singletons total:   {sa_sing}")
    lines.append("")
    lines.append(f"Possible mislabelled members (SA taxon, international gene ID): {len(mismatches)}")
    lines.append(f"   --fix-labels was {'ON' if args.fix_labels else 'OFF'}")
    if bad_lines:
        lines.append(f"Skipped malformed group lines: {bad_lines}")

    report = "\n".join(lines)
    print(report)
    with open(os.path.join(args.outdir, "summary.txt"), "w") as o:
        o.write(report + "\n")
    with open(os.path.join(args.outdir, "mislabelled_members.txt"), "w") as o:
        for gid, m in mismatches:
            o.write(f"{gid}\t{m}\n")


if __name__ == "__main__":
    sys.exit(main())
