"""
make_wc_triplets_ndx.py

Builds a GROMACS index file with Watson-Crick H-bond triplet groups for
circular DNA, in H-D-A order (hydrogen, donor heavy atom, acceptor),
for use with gmx angle.

Assumes the standard circular numbering: strand 1 is residues 1..N_bp,
strand 2 is (N_bp+1)..(2*N_bp), and residue i pairs with residue
(2*N_bp + 1 - i).

Donors: DT N3-H3 -> DA N1, and DG N1-H1 -> DC N3

Usage for an 80 bp circular DNA (160 residues total):
    python3 make_wc_triplets_ndx.py structure.gro 80 [output_dir]
"""

import sys
from pathlib import Path

bases = {"DA", "DT", "DG", "DC"}


def core(resname):
    # collapse terminal labels like DA3/DT5 down to the plain base name
    rn = resname.strip().upper()
    if rn in bases:
        return rn
    if rn[-1:] in {"3", "5"} and rn[:-1] in bases:
        return rn[:-1]
    return rn


def parse_gro(gro):
    lines = gro.read_text().splitlines()
    natoms = int(lines[1].strip())
    lines = lines[2:2 + natoms]
    atoms = []
    for ln in lines:
        resid = int(ln[0:5])
        resname = ln[5:10].strip()
        aname = ln[10:15].strip()
        idx = int(ln[15:20])
        if core(resname) in bases:
            atoms.append((resid, resname, aname, idx))
    if not atoms:
        raise ValueError("No DNA bases (DA/DT/DG/DC) found in GRO")
    return atoms


def lookup_table(atoms):
    # (resid, base, atom_name) -> atom index
    table = {}
    for resid, resname, aname, idx in atoms:
        b = core(resname)
        table[(resid, b, aname.strip().upper())] = idx
    return table


def get_base_type(atoms, resid):
    for r, resname, _, _ in atoms:
        if r == resid:
            return core(resname)
    return None


def triplet_for_pair(table, res_a, base_a, res_b, base_b):
    # returns H-D-A atom indices, or None if anything's missing/non-WC
    if base_a == "DA" and base_b == "DT":
        # donor is DT on strand b
        h = table.get((res_b, "DT", "H3"))
        d = table.get((res_b, "DT", "N3"))
        a = table.get((res_a, "DA", "N1"))
    elif base_a == "DT" and base_b == "DA":
        # donor is DT on strand a
        h = table.get((res_a, "DT", "H3"))
        d = table.get((res_a, "DT", "N3"))
        a = table.get((res_b, "DA", "N1"))
    elif base_a == "DG" and base_b == "DC":
        # donor is DG on strand a
        h = table.get((res_a, "DG", "H1"))
        d = table.get((res_a, "DG", "N1"))
        a = table.get((res_b, "DC", "N3"))
    elif base_a == "DC" and base_b == "DG":
        # donor is DG on strand b
        h = table.get((res_b, "DG", "H1"))
        d = table.get((res_b, "DG", "N1"))
        a = table.get((res_a, "DC", "N3"))
    else:
        return None

    if None in (h, d, a):
        return None
    return [h, d, a]


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 make_wc_triplets_ndx_circular.py structure.gro N_bp [output_dir]")
        print("  N_bp = number of base pairs (e.g. 80 for 160-residue circular DNA)")
        print("  output_dir = optional, defaults to the input gro's directory")
        sys.exit(1)

    gro_path = Path(sys.argv[1])
    n_bp = int(sys.argv[2])

    if len(sys.argv) >= 4:
        output_dir = Path(sys.argv[3])
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = gro_path.parent

    out_path = output_dir / "wc_triplets.ndx"

    atoms = parse_gro(gro_path)
    table = lookup_table(atoms)

    residues = sorted(set(r for r, _, _, _ in atoms))
    print(f"Found {len(residues)} residues: {min(residues)} to {max(residues)}")

    groups = []
    all_triplets = []

    # chain 1 residue i pairs with chain 2 residue (2*n_bp + 1 - i)
    for i in range(1, n_bp + 1):
        res_a = i
        res_b = (2 * n_bp + 1) - i

        base_a = get_base_type(atoms, res_a)
        base_b = get_base_type(atoms, res_b)

        if not base_a or not base_b:
            print(f"warning: missing base for pair {i}: res {res_a} ({base_a}) <-> res {res_b} ({base_b})")
            continue

        if base_a not in ("DA", "DT", "DG", "DC") or base_b not in ("DA", "DT", "DG", "DC"):
            print(f"warning: non-WC pair {i}: {base_a} (res {res_a}) <-> {base_b} (res {res_b})")
            continue

        triplet = triplet_for_pair(table, res_a, base_a, res_b, base_b)

        if triplet:
            groups.append((f"BP{i:02d}_trip", triplet))
            all_triplets.extend(triplet)
            print(f"triplet {i:02d}: res{res_a:3d} ({base_a}) <-> res{res_b:3d} ({base_b}) | atoms {triplet[0]:5d} - {triplet[1]:5d} - {triplet[2]:5d}")
        else:
            print(f"warning: missing atoms or non-WC pair {i}: res {res_a} ({base_a}) <-> res {res_b} ({base_b})")

    with out_path.open("w") as f:
        for name, idxs in groups:
            f.write(f"[ {name} ]\n")
            f.write(" ".join(map(str, idxs)) + "\n\n")

        if all_triplets:
            f.write("[ all_triplets ]\n")
            for k in range(0, len(all_triplets), 3):
                f.write(f"{all_triplets[k]:>5} {all_triplets[k+1]:>5} {all_triplets[k+2]:>5}\n")
            f.write("\n")

    print(f"\nwrote {len(groups)} triplet groups and [ all_triplets ] to {out_path}")


if __name__ == "__main__":
    main()    for prefix in ("DA", "DT", "DG", "DC", "ADE", "THY", "GUA", "CYT", "A", "T", "G", "C"):
        if letters.startswith(prefix):
            return prefix
    return letters


def base_letter_from_resname(rn_norm):
    if rn_norm.startswith("DA") or rn_norm in ("ADE", "A"):
        return "A"
    if rn_norm.startswith("DT") or rn_norm in ("THY", "T"):
        return "T"
    if rn_norm.startswith("DG") or rn_norm in ("GUA", "G"):
        return "G"
    if rn_norm.startswith("DC") or rn_norm in ("CYT", "C"):
        return "C"
    return None


def parse_gro_atoms(gro_path):
    # returns res_order (resid, resname) in file order, deduped per residue,
    # and idx_by_resid_atom mapping (resid, atomname) -> 1-based atom index
    with open(gro_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    if len(lines) < 3:
        raise ValueError("GRO file looks truncated")

    natoms = int(lines[1].strip())
    atom_lines = lines[2:2 + natoms]

    res_order = []
    seen_resids = set()
    idx_by_resid_atom = {}

    for i, line in enumerate(atom_lines, start=1):
        resid = int(line[0:5].strip())
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        if resid not in seen_resids:
            res_order.append((resid, resname))
            seen_resids.add(resid)
        idx_by_resid_atom[(resid, atomname)] = i

    return res_order, idx_by_resid_atom


def build_dna_residue_list(res_order):
    dna_resids, dna_bases = [], []
    for resid, rn in res_order:
        base = base_letter_from_resname(normalize_resname(rn))
        if base is not None:
            dna_resids.append(resid)
            dna_bases.append(base)
    if not dna_resids:
        raise ValueError("No DA/DT/DG/DC residues found in the gro")
    return dna_resids, dna_bases


def lookup_idx(idx_by_resid_atom, resid, aname):
    # exact name first, then case variants, then a slow scan for
    # whitespace weirdness in the gro
    if (resid, aname) in idx_by_resid_atom:
        return idx_by_resid_atom[(resid, aname)]
    if (resid, aname.upper()) in idx_by_resid_atom:
        return idx_by_resid_atom[(resid, aname.upper())]
    if (resid, aname.lower()) in idx_by_resid_atom:
        return idx_by_resid_atom[(resid, aname.lower())]
    for (r, a), idx in idx_by_resid_atom.items():
        if r == resid and a.strip().upper() == aname.upper():
            return idx
    raise KeyError(f"couldn't find atom '{aname}' on resid {resid}")


def triplet_for_pair_dnaorder(o1, o2, dna_resids, dna_bases, idx_by_resid_atom):
    # o1, o2 are 1-based positions in the DNA-only list
    n = len(dna_resids)
    if not (1 <= o1 <= n and 1 <= o2 <= n):
        raise ValueError(f"pair ({o1},{o2}) is outside DNA range 1..{n}")

    r1, b1 = dna_resids[o1 - 1], dna_bases[o1 - 1]
    r2, b2 = dna_resids[o2 - 1], dna_bases[o2 - 1]

    donors, acceptors = {"G", "T"}, {"A", "C"}
    if b1 in donors and b2 in acceptors:
        donor_res, donor_base, accept_res, accept_base = r1, b1, r2, b2
    elif b2 in donors and b1 in acceptors:
        donor_res, donor_base, accept_res, accept_base = r2, b2, r1, b1
    else:
        raise ValueError(f"pair {o1}({b1})-{o2}({b2}) isn't a valid donor/acceptor combo")

    if donor_base == "G":
        if accept_base != "C":
            raise ValueError(f"G donor needs a C acceptor, got {accept_base}")
        H = lookup_idx(idx_by_resid_atom, donor_res, "H1")
        D = lookup_idx(idx_by_resid_atom, donor_res, "N1")
        A = lookup_idx(idx_by_resid_atom, accept_res, "N3")
        return (H, D, A)

    # donor_base == "T"
    if accept_base != "A":
        raise ValueError(f"T donor needs an A acceptor, got {accept_base}")
    H = lookup_idx(idx_by_resid_atom, donor_res, "H3")
    D = lookup_idx(idx_by_resid_atom, donor_res, "N3")
    A = lookup_idx(idx_by_resid_atom, accept_res, "N1")
    return (H, D, A)


def read_pairs_dnaorder(pairs_path):
    pairs = []
    with open(pairs_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) >= 2:
                pairs.append((int(parts[0]), int(parts[1])))
    if not pairs:
        raise ValueError("no pairs found in pairs.tsv")
    return pairs


def write_ndx(out_path, triplets):
    with open(out_path, "w", encoding="utf-8") as f:
        for k, (h, d, a) in enumerate(triplets, start=1):
            f.write(f"[ BP{str(k).zfill(2)} ]\n")
            f.write(f"{h:6d} {d:6d} {a:6d}\n\n")


def main():
    if len(sys.argv) != 4:
        print("Usage: python make_wc_triplets_ndx.py structure.gro pairs.tsv wc_triplets.ndx")
        sys.exit(1)

    gro_path = Path(sys.argv[1])
    pairs_path = Path(sys.argv[2])
    out_path = Path(sys.argv[3])

    res_order, idx_by_resid_atom = parse_gro_atoms(gro_path)
    dna_resids, dna_bases = build_dna_residue_list(res_order)
    pairs = read_pairs_dnaorder(pairs_path)

    max_id = max(max(a, b) for a, b in pairs)
    if max_id > len(dna_resids):
        raise ValueError(f"pairs.tsv references DNA index {max_id} but only {len(dna_resids)} DNA residues found")

    triplets = [triplet_for_pair_dnaorder(o1, o2, dna_resids, dna_bases, idx_by_resid_atom) for o1, o2 in pairs]
    write_ndx(out_path, triplets)
    print(f"wrote {len(triplets)} groups to {out_path}")


if __name__ == "__main__":
    main()
