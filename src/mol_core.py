"""
Core parsing and fragment-mass utilities for the rearrangement classifier.
Pure standard library -- no RDKit/OpenBabel dependency (none available offline).
"""

import re
from collections import defaultdict, deque

# ---------------------------------------------------------------------------
# Monoisotopic atomic masses (Da) for elements expected in METLIN structures
# ---------------------------------------------------------------------------
MONO_MASS = {
    'H': 1.0078250319, 'C': 12.0000000, 'N': 14.0030740052,
    'O': 15.9949146221, 'F': 18.9984032, 'Na': 22.98976928,
    'Si': 27.9769265327, 'P': 30.97376151, 'S': 31.97207069,
    'Cl': 34.96885268, 'K': 38.9637069, 'Br': 78.9183371,
    'I': 126.904473,
}
PROTON_MASS = 1.00727646688
ELECTRON_MASS = 0.00054858

# Standard valence used for implicit-H bookkeeping (neutral, common oxidation state)
VALENCE = {
    'C': 4, 'N': 3, 'O': 2, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
    'S': 2, 'P': 3, 'Si': 4, 'Na': 1, 'K': 1, 'H': 1,
}


class Atom:
    __slots__ = ('idx', 'element', 'x', 'y', 'z')

    def __init__(self, idx, element, x, y, z):
        self.idx = idx  # 0-indexed
        self.element = element
        self.x, self.y, self.z = x, y, z


class Bond:
    __slots__ = ('a1', 'a2', 'order')

    def __init__(self, a1, a2, order):
        self.a1 = a1  # 0-indexed
        self.a2 = a2
        self.order = order  # 1,2,3 int; 4 = aromatic (treated as 1.5 for valence math)


class Molecule:
    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds
        self.adj = defaultdict(list)  # atom_idx -> list of (neighbor_idx, bond_idx)
        for bi, b in enumerate(bonds):
            self.adj[b.a1].append((b.a2, bi))
            self.adj[b.a2].append((b.a1, bi))

    def bond_order_sum(self, atom_idx, exclude_bond_idx=None, exclude_bond_idxs=None):
        exclude_set = set()
        if exclude_bond_idx is not None:
            exclude_set.add(exclude_bond_idx)
        if exclude_bond_idxs is not None:
            exclude_set.update(exclude_bond_idxs)
        total = 0.0
        for (_, bi) in self.adj[atom_idx]:
            if bi in exclude_set:
                continue
            order = self.bonds[bi].order
            total += 1.5 if order == 4 else order
        return total

    def is_bridge(self, bond_idx):
        """A bond is a bridge (non-ring) if removing it disconnects its two endpoints."""
        b = self.bonds[bond_idx]
        seen = {b.a1}
        dq = deque([b.a1])
        while dq:
            cur = dq.popleft()
            for (nbr, bi) in self.adj[cur]:
                if bi == bond_idx:
                    continue
                if nbr not in seen:
                    seen.add(nbr)
                    dq.append(nbr)
        return b.a2 not in seen

    def component_after_removal(self, bond_idx, start_atom):
        """BFS to find all atoms reachable from start_atom with bond_idx removed."""
        seen = {start_atom}
        dq = deque([start_atom])
        while dq:
            cur = dq.popleft()
            for (nbr, bi) in self.adj[cur]:
                if bi == bond_idx:
                    continue
                if nbr not in seen:
                    seen.add(nbr)
                    dq.append(nbr)
        return seen

    def components_after_removal(self, bond_idxs):
        """Remove a SET of bonds (any size) and return ALL resulting connected
        components (list of atom-index sets) across the whole molecule --
        not just the component containing one starting atom. May return more
        than 2 components (caller decides what to do with that)."""
        exclude = set(bond_idxs)
        all_atoms = set(range(len(self.atoms)))
        unvisited = set(all_atoms)
        components = []
        while unvisited:
            start = next(iter(unvisited))
            seen = {start}
            dq = deque([start])
            while dq:
                cur = dq.popleft()
                for (nbr, bi) in self.adj[cur]:
                    if bi in exclude:
                        continue
                    if nbr not in seen:
                        seen.add(nbr)
                        dq.append(nbr)
            components.append(seen)
            unvisited -= seen
        return components

    def ring_systems(self):
        """Identify fused-ring clusters: connected groups of non-bridge (ring)
        bonds, where two ring bonds belong to the same system if they share
        an atom (directly, or via a chain of other ring bonds). A single
        isolated ring is one system; a fused bicyclic/polycyclic structure
        is also one system, since its rings share bonds/atoms transitively.
        Returns a list of bond-index sets, one per distinct ring system."""
        ring_bond_idxs = [bi for bi in range(len(self.bonds)) if not self.is_bridge(bi)]
        ring_bond_set = set(ring_bond_idxs)
        # Build bond-level adjacency: two ring bonds are linked if they share an atom
        atom_to_ring_bonds = defaultdict(list)
        for bi in ring_bond_idxs:
            b = self.bonds[bi]
            atom_to_ring_bonds[b.a1].append(bi)
            atom_to_ring_bonds[b.a2].append(bi)

        visited = set()
        systems = []
        for bi in ring_bond_idxs:
            if bi in visited:
                continue
            seen = {bi}
            dq = deque([bi])
            visited.add(bi)
            while dq:
                cur = dq.popleft()
                b = self.bonds[cur]
                for atom in (b.a1, b.a2):
                    for other_bi in atom_to_ring_bonds[atom]:
                        if other_bi not in seen:
                            seen.add(other_bi)
                            visited.add(other_bi)
                            dq.append(other_bi)
            systems.append(seen)
        return systems


def parse_molblock(lines):
    """
    Parse a V2000 MOL block (list of lines, no leading '$$$$' or trailing blank).
    lines[0:3] = header/title/comment; lines[3] = counts line.
    Robust to counts-line spacing variations across ChemDraw/Marvin/RDKit writers.
    """
    counts_line = lines[3]
    # Counts line is fixed-width columns (3 chars each for natoms/nbonds), NOT
    # whitespace-delimited -- e.g. "93100  0  0..." is 93 atoms + 100 bonds with
    # no space between them. A plain .split() misparses this whenever a count
    # is exactly 3 digits and butts up against its neighbor.
    natoms = int(counts_line[0:3])
    nbonds = int(counts_line[3:6])

    atoms = []
    for i in range(natoms):
        line = lines[4 + i]
        x = float(line[0:10])
        y = float(line[10:20])
        z = float(line[20:30])
        element = line[31:34].strip()
        atoms.append(Atom(i, element, x, y, z))

    bonds = []
    bond_start = 4 + natoms
    for i in range(nbonds):
        line = lines[bond_start + i]
        a1 = int(line[0:3]) - 1  # MOL is 1-indexed
        a2 = int(line[3:6]) - 1
        order = int(line[6:9])
        bonds.append(Bond(a1, a2, order))

    return Molecule(atoms, bonds)


def split_sdf_records(path):
    """Yield (molblock_lines, properties_dict) for each record in an SDF file."""
    with open(path, 'r', errors='replace') as f:
        buf = []
        for line in f:
            if line.rstrip('\n') == '$$$$':
                yield _parse_record(buf)
                buf = []
            else:
                buf.append(line.rstrip('\n'))
        if buf:
            yield _parse_record(buf)


def _parse_record(buf):
    # Find where the MOL block ends (the "M  END" line)
    end_idx = None
    for i, line in enumerate(buf):
        # Some records have a stray trailing '"' appended directly to the
        # M END line with no separating whitespace (e.g. 'M  END"') -- an
        # artifact of the multi-source SDF assembly, not a real structural
        # difference. Strip trailing quote chars before the exact-match check.
        if line.strip().rstrip('"') == 'M  END':
            end_idx = i
            break
    if end_idx is None:
        return None, {}
    molblock_lines = buf[:end_idx + 1]

    props = {}
    rest = buf[end_idx + 1:]
    field_re = re.compile(r'^>\s*<(.+?)>')
    i = 0
    while i < len(rest):
        m = field_re.match(rest[i])
        if m:
            key = m.group(1)
            i += 1
            val_lines = []
            while i < len(rest) and rest[i].strip() != '':
                val_lines.append(rest[i])
                i += 1
            props[key] = '\n'.join(val_lines)
        else:
            i += 1
    return molblock_lines, props


def formula_and_mass(element_counts):
    """Given a dict {element: count}, return (formula_str, monoisotopic_mass)."""
    mass = 0.0
    for el, n in element_counts.items():
        if n <= 0:
            continue
        if el not in MONO_MASS:
            return None, None  # unsupported element, skip
        mass += MONO_MASS[el] * n
    # Hill system formula string (C first, H second, rest alphabetical)
    parts = []
    if 'C' in element_counts and element_counts['C'] > 0:
        parts.append(f"C{element_counts['C'] if element_counts['C'] != 1 else ''}")
    if 'H' in element_counts and element_counts['H'] > 0:
        parts.append(f"H{element_counts['H'] if element_counts['H'] != 1 else ''}")
    for el in sorted(k for k in element_counts if k not in ('C', 'H')):
        n = element_counts[el]
        if n > 0:
            parts.append(f"{el}{n if n != 1 else ''}")
    return ''.join(parts), mass


def fragment_element_counts(mol, atom_set, removed_bond_idx, cleave_atom):
    """
    Element counts (including implicit H) for the fragment consisting of atom_set,
    where removed_bond_idx was the bond cut, and cleave_atom is the atom in this
    fragment that was directly bonded across the break (gets the valence cap).
    Returns dict element->count for ΔH=0 baseline (simple homolytic cap).
    """
    counts = defaultdict(int)
    removed_order_at_cleave = None
    b = mol.bonds[removed_bond_idx]
    removed_order = 1.5 if b.order == 4 else b.order

    for aidx in atom_set:
        atom = mol.atoms[aidx]
        el = atom.element
        if el not in VALENCE:
            continue
        counts[el] += 1
        bond_sum = mol.bond_order_sum(aidx, exclude_bond_idx=None)
        # For the cleave atom, the removed bond is not in adj anymore conceptually,
        # but bond_order_sum already only sums bonds still present in mol.bonds;
        # the removed bond IS still in mol.bonds (we didn't delete it), so exclude it:
        if aidx == cleave_atom:
            bond_sum = mol.bond_order_sum(aidx, exclude_bond_idx=removed_bond_idx)
            implicit_h = VALENCE[el] - bond_sum  # true zero-rearrangement anchor, no extra cap
        else:
            implicit_h = VALENCE[el] - bond_sum
        implicit_h = max(0, round(implicit_h))
        counts['H'] += implicit_h

    return dict(counts)
