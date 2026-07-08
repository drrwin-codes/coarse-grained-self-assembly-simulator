"""
Identify amphiphile aggregates and classify their shape from the gyration
tensor. An aggregate is defined by chain bonds + tail-tail lattice contacts
ONLY. Head-head or head-tail (cross-chain) lattice adjacency does NOT link
two aggregates -- otherwise crowded head groups from separate micelles that
merely touch get counted as one giant merged cluster, which is a labeling
artifact, not real coalescence. This matches how aggregation number /
core connectivity is defined in the actual surfactant self-assembly
literature (aggregates = shared hydrophobic core).
"""
import numpy as np


def _min_image_delta(x1, y1, x2, y2, L):
    dx, dy = x2 - x1, y2 - y1
    if dx > L / 2: dx -= L
    if dx < -L / 2: dx += L
    if dy > L / 2: dy -= L
    if dy < -L / 2: dy += L
    return dx, dy


def find_aggregates(system):
    L = system.L
    grid = system.grid
    chains = system.chains
    type_seq = system.type_seq

    # site -> (chain_idx, bead_idx)
    site_to_bead = {}
    for ci, chain in enumerate(chains):
        for bi, site in enumerate(chain):
            site_to_bead[site] = (ci, bi)

    deltas = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    visited = set()
    aggregates = []

    for ci, chain in enumerate(chains):
        for bi in range(len(chain)):
            if (ci, bi) in visited:
                continue
            # BFS over bead-graph: bond edges (always) + tail-tail contact edges
            stack = [(ci, bi, float(chain[bi][0]), float(chain[bi][1]))]
            visited.add((ci, bi))
            members = []  # (chain_idx, bead_idx, ux, uy)
            while stack:
                cci, bbi, ux, uy = stack.pop()
                members.append((cci, bbi, ux, uy))
                cur_site = chains[cci][bbi]
                cur_type = type_seq[bbi]

                # bond neighbors within the same chain
                for nb_bi in (bbi - 1, bbi + 1):
                    if 0 <= nb_bi < len(chains[cci]) and (cci, nb_bi) not in visited:
                        nb_site = chains[cci][nb_bi]
                        dx, dy = _min_image_delta(cur_site[0], cur_site[1],
                                                   nb_site[0], nb_site[1], L)
                        visited.add((cci, nb_bi))
                        stack.append((cci, nb_bi, ux + dx, uy + dy))

                # tail-tail contact neighbors (cross-chain allowed), matching
                # the diagonal-inclusive interaction range used in energetics
                if cur_type == 'T':
                    contact_deltas = [(1, 0), (-1, 0), (0, 1), (0, -1),
                                       (1, 1), (1, -1), (-1, 1), (-1, -1)]
                    for dx0, dy0 in contact_deltas:
                        nx, ny = (cur_site[0] + dx0) % L, (cur_site[1] + dy0) % L
                        if grid[nx, ny] != 'T':
                            continue
                        nb = site_to_bead.get((nx, ny))
                        if nb is None or nb in visited:
                            continue
                        visited.add(nb)
                        stack.append((nb[0], nb[1], ux + dx0, uy + dy0))
            aggregates.append(members)
    return aggregates


def analyze_aggregate(members, L):
    pts = np.array([(m[2], m[3]) for m in members])
    chain_ids = set(m[0] for m in members)
    n_beads = len(pts)
    agg_number = len(chain_ids)  # number of distinct chains in this aggregate

    com = pts.mean(axis=0)
    rel = pts - com
    S = (rel.T @ rel) / n_beads
    eigvals = np.linalg.eigvalsh(S)
    lam2, lam1 = eigvals[0], eigvals[1]
    ecc = 0.0 if lam1 <= 1e-9 else 1.0 - (lam2 / lam1)

    span_x = pts[:, 0].max() - pts[:, 0].min()
    span_y = pts[:, 1].max() - pts[:, 1].min()
    spans = (span_x >= L - 1) or (span_y >= L - 1)

    if spans:
        label = 'bilayer/vesicle'
    elif agg_number <= 6:
        # at this size, single-snapshot eccentricity is dominated by
        # thermal noise (empirically spans ~0.06-0.98 regardless of shape) --
        # aggregation number alone is the stable signal here.
        label = 'micelle'
    elif agg_number <= 10:
        label = 'micelle' if ecc < 0.55 else 'cylindrical/wormlike'
    else:
        label = 'cylindrical/wormlike'

    return {
        'n_beads': n_beads, 'agg_number': agg_number,
        'eccentricity': ecc, 'spans_box': bool(spans), 'label': label,
    }


def find_water_pockets(system):
    """Flood-fill the solvent (periodic BC). The largest water cluster is
    'bulk' solvent; any other, smaller, disconnected water cluster is a
    pocket trapped inside a closed amphiphile ring -- i.e. a real vesicle
    interior, not just an elongated/high-eccentricity aggregate."""
    L = system.L
    grid = system.grid
    visited = set()
    deltas = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    clusters = []
    for x0 in range(L):
        for y0 in range(L):
            if grid[x0, y0] != 'W' or (x0, y0) in visited:
                continue
            stack = [(x0, y0)]
            visited.add((x0, y0))
            sites = []
            while stack:
                x, y = stack.pop()
                sites.append((x, y))
                for dx, dy in deltas:
                    nx, ny = (x + dx) % L, (y + dy) % L
                    if grid[nx, ny] == 'W' and (nx, ny) not in visited:
                        visited.add((nx, ny))
                        stack.append((nx, ny))
            clusters.append(sites)
    clusters.sort(key=len, reverse=True)
    bulk = clusters[0] if clusters else []
    pockets = clusters[1:]
    return bulk, pockets


def find_enclosing_aggregate(pocket_sites, system, aggregates):
    """Which aggregate encloses this water pocket? Returns the aggregate
    index only if the pocket is bordered EXCLUSIVELY by that one aggregate --
    a genuine closed ring around its own interior. If more than one distinct
    aggregate borders the pocket, it's a crowding gap squeezed between
    separate structures, not a vesicle, and this returns None."""
    L = system.L
    site_to_agg = {}
    for idx, members in enumerate(aggregates):
        for (ci, bi, ux, uy) in members:
            site_to_agg[system.chains[ci][bi]] = idx
    deltas = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    bordering = set()
    for (x, y) in pocket_sites:
        for dx, dy in deltas:
            nx, ny = (x + dx) % L, (y + dy) % L
            agg_idx = site_to_agg.get((nx, ny))
            if agg_idx is not None:
                bordering.add(agg_idx)
    if len(bordering) != 1:
        return None
    return next(iter(bordering))


MIN_POCKET_SIZE = 4   # rules out single/double-site lattice-roughness notches
MIN_AGG_FOR_VESICLE = 10  # a handful of chains can't physically close a ring


def classify_system_with_vesicles(system):
    """Full classification pass: shape-based label, then override with
    'vesicle' only for aggregates that (a) exclusively enclose a water
    pocket, (b) the pocket is big enough to be a real interior rather than
    lattice-scale surface roughness, and (c) the enclosing aggregate has
    enough chains to plausibly form a closed ring at all."""
    aggregates = find_aggregates(system)
    results = [analyze_aggregate(m, system.L) for m in aggregates]
    bulk, pockets = find_water_pockets(system)
    vesicle_info = []
    for pocket in pockets:
        if len(pocket) < MIN_POCKET_SIZE:
            continue
        agg_idx = find_enclosing_aggregate(pocket, system, aggregates)
        if agg_idx is None:
            continue
        if results[agg_idx]['agg_number'] < MIN_AGG_FOR_VESICLE:
            continue
        results[agg_idx]['label'] = 'vesicle'
        results[agg_idx]['pocket_size'] = len(pocket)
        vesicle_info.append({'agg_idx': agg_idx, 'pocket_size': len(pocket)})
    return results, vesicle_info



    aggregates = find_aggregates(system)
    return [analyze_aggregate(m, system.L) for m in aggregates]


def summarize(results):
    from collections import Counter
    counts = Counter(r['label'] for r in results)
    agg_numbers = [r['agg_number'] for r in results]
    return counts, agg_numbers
