
import argparse
import os
from collections import Counter

DEFAULT_GROUPS = "C:/HonoursProj/ortholog_clustering/groups.txt"
DEFAULT_SINGLETONS = "C:/HonoursProj/ortholog_clustering/singletons.txt"

SA_VARIETIES = {
    "erica": "Erica",
    "kadie": "Kadie",
    "malgas": "Malgas",
    "nemesia": "Nemesia",
    "rgt_planet": "RGT_Planet",
    "rgtplanet": "RGT_Planet",
    "planet": "RGT_Planet",
    "uvc8": "UVC8"
}

SA5 = {"Erica", "Kadie", "Malgas", "Nemesia", "RGT_Planet"}
SA6 = SA5 | {"UVC8"}

INTL_MARKERS = ["MOREX", "IGRI", "HOCKETT", "GOLDEN_PROMISE", "BARKE", "MAXIMUS"]


def normalise(variety):
    return SA_VARIETIES.get(variety.lower(), variety)


def parse_variety(member, fix_labels, mismatches, gid):
    variety, _, gene = member.partition("|")
    variety = normalise(variety)

    if variety in SA6 and any(marker in gene.upper() for marker in INTL_MARKERS):
        mismatches.append((gid, member))

        if fix_labels:
            variety = "INTERNATIONAL_MISLABELLED"

    return variety


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("groups", nargs="?", default=DEFAULT_GROUPS)
    parser.add_argument("singletons", nargs="?", default=DEFAULT_SINGLETONS)
    parser.add_argument("--outdir", default="orthomcl_results")
    parser.add_argument("--fix-labels", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    categories = {
        "sa_only_all6": [],
        "sa5_only_no_uvc8": [],
        "rgt_no_international": [],
        "rgt_planet_only": [],
        "uvc8_only": []
    }

    total_groups = 0
    variety_counts = Counter()
    mismatches = []
    group_lines = {}
    bad_lines = 0

    # Read OrthoMCL groups
    with open(args.groups) as f, open(os.path.join(args.outdir, "group_taxa_matrix.tsv"), "w") as matrix:
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

            varieties = {
                parse_variety(m, args.fix_labels, mismatches, gid)
                for m in members
            }

            for variety in varieties:
                variety_counts[variety] += 1

            matrix.write(f"{gid}\t{len(members)}\t{','.join(sorted(varieties))}\n")

            # Identify South African groups
            if varieties <= SA6:
                group_lines[gid] = " ".join(members)
                categories["sa_only_all6"].append(gid)

                if "UVC8" not in varieties:
                    categories["sa5_only_no_uvc8"].append(gid)

                if "RGT_Planet" in varieties:
                    categories["rgt_no_international"].append(gid)

                if varieties == {"RGT_Planet"}:
                    categories["rgt_planet_only"].append(gid)

                if varieties == {"UVC8"}:
                    categories["uvc8_only"].append(gid)

    # Save category files
    for name, groups in categories.items():
        with open(os.path.join(args.outdir, name + ".txt"), "w") as f:
            for gid in groups:
                f.write(f"{gid}: {group_lines[gid]}\n")

    # Count singletons per variety
    singleton_counts = Counter()

    with open(args.singletons) as f:
        for line in f:
            for item in line.split():
                variety, _, gene = item.partition("|")
                variety = normalise(variety)

                if variety in SA6 and any(m in gene.upper() for m in INTL_MARKERS):
                    if args.fix_labels:
                        variety = "INTERNATIONAL_MISLABELLED"

                singleton_counts[variety] += 1

    # Generate summary
    report = [
        f"Total groups:                                   {total_groups}",
        f"Total singletons:                               {sum(singleton_counts.values())}",
        "",
        f"SA-only groups (all 6 incl. UVC8):              {len(categories['sa_only_all6'])}",
        f"Unique to 5 SA varieties (no UVC8):             {len(categories['sa5_only_no_uvc8'])}",
        f"Contain RGT_Planet, no international:           {len(categories['rgt_no_international'])}",
        f"   of which RGT_Planet only:                    {len(categories['rgt_planet_only'])}",
        f"UVC8-only groups:                               {len(categories['uvc8_only'])}",
        "",
        "Groups containing each variety:"
    ]

    for variety, count in variety_counts.most_common():
        report.append(f"   {variety:28s}{count}")

    report.append("")
    report.append("Singletons per variety:")

    for variety, count in singleton_counts.most_common():
        report.append(f"   {variety:28s}{count}")

    sa_singletons = sum(singleton_counts[v] for v in SA6)
    report.append(f"   SA-grown singletons total:   {sa_singletons}")

    report.append("")
    report.append(f"Possible mislabelled members: {len(mismatches)}")
    report.append(f"--fix-labels was {'ON' if args.fix_labels else 'OFF'}")

    if bad_lines:
        report.append(f"Skipped malformed group lines: {bad_lines}")

    report = "\n".join(report)
    print(report)

    with open(os.path.join(args.outdir, "summary.txt"), "w") as f:
        f.write(report + "\n")

    with open(os.path.join(args.outdir, "mislabelled_members.txt"), "w") as f:
        for gid, member in mismatches:
            f.write(f"{gid}\t{member}\n")


if __name__ == "__main__":
    main()

