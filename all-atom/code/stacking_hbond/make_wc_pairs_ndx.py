"""
make_wc_pairs_ndx.py

Builds a GROMACS index file with Watson-Crick H-bond pairs for circular
DNA. Assumes the standard circular numbering where chain 1 (residues
1..N_bp) pairs with chain 2 in reverse: residue i pairs with
2*N_bp + 1 - i.

Usage:
    python3 make_wc_pairs_ndx.py structure.gro 80 [output_dir]

80 here is N_bp, the number of base pairs (so total DNA residues = 2*N_bp).
output_dir is optional, defaults to wherever the gro file lives.
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


def write_group(f, name, values, per_line=15):
    f.write(f"[ {name} ]\n")
    for idx, atom in enumerate(values):
        f.write(f"{atom:6d}")
        f.write("\n" if (idx + 1) % per_line == 0 else " ")
    if len(values) % per_line != 0:
        f.write("\n")
    f.write("\n")


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 make_wc_pairs_ndx.py structure.gro N_bp [output_dir]")
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

    pairs_out = output_dir / "wc_pairs.ndx"

    atoms = parse_gro(gro_path)
    table = lookup_table(atoms)

    residues = sorted(set(r for r, _, _, _ in atoms))
    print(f"Found {len(residues)} residues: {min(residues)} to {max(residues)}")

    pair_groups = []
    all_pairs = []
    donors = []
    acceptors = []

    # chain 1 residue i pairs with chain 2 residue (2*n_bp + 1 - i)
    for i in range(1, n_bp + 1):
        res_a = i
        res_b = (2 * n_bp + 1) - i

        base_a = get_base_type(atoms, res_a)
        base_b = get_base_type(atoms, res_b)

        if not base_a or not base_b:
            print(f"warning: missing base for pair {i}: res {res_a} ({base_a}) <-> res {res_b} ({base_b})")
            continue

        donor, acceptor = None, None

        if base_a == "DA" and base_b == "DT":
            acceptor = table.get((res_a, "DA", "N1"))
            donor = table.get((res_b, "DT", "N3"))
        elif base_a == "DT" and base_b == "DA":
            donor = table.get((res_a, "DT", "N3"))
            acceptor = table.get((res_b, "DA", "N1"))
        elif base_a == "DG" and base_b == "DC":
            donor = table.get((res_a, "DG", "N1"))
            acceptor = table.get((res_b, "DC", "N3"))
        elif base_a == "DC" and base_b == "DG":
            acceptor = table.get((res_a, "DC", "N3"))
            donor = table.get((res_b, "DG", "N1"))
        else:
            print(f"warning: non-WC pair {i}: {base_a} (res {res_a}) <-> {base_b} (res {res_b})")
            continue

        if donor and acceptor:
            pair = [donor, acceptor]
            name = f"P{i:02d}"
            pair_groups.append((name, pair))
            all_pairs.extend(pair)
            donors.append(donor)
            acceptors.append(acceptor)
            print(f"pair {i:02d}: res{res_a:3d} ({base_a}) <-> res{res_b:3d} ({base_b}) | atoms {donor:5d} - {acceptor:5d}")
        else:
            print(f"warning: missing atoms for pair {i}: res {res_a} ({base_a}) <-> res {res_b} ({base_b})")

    with pairs_out.open("w") as f:
        for name, pair in pair_groups:
            f.write(f"[ {name} ]\n{pair[0]:6d} {pair[1]:6d}\n\n")

        if all_pairs:
            write_group(f, "all_pairs", all_pairs)
        if donors:
            write_group(f, "donors", donors)
        if acceptors:
            write_group(f, "acceptors", acceptors)

    print(f"\nwrote {len(pair_groups)} pairs to {pairs_out}")
    print(f"  all_pairs: {len(all_pairs)} atoms")
    print(f"  donors: {len(donors)} atoms")
    print(f"  acceptors: {len(acceptors)} atoms")


if __name__ == "__main__":
    main()
