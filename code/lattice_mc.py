"""
Minimal 2D lattice Monte Carlo model of amphiphile self-assembly.

Amphiphiles are block chains of H (head) and T (tail) beads on a
periodic square lattice. Water is implicit (empty sites). Interactions
are nearest-neighbor only. No shape is imposed -- whatever aggregate
geometry forms, forms on its own from the energy rules + Metropolis MC.
"""
import numpy as np
import random

# --- interaction parameters (kT = 1 units) ---
EPS = {
    ('T', 'T'): -1.8, ('W', 'W'): 0.0, ('H', 'H'): 0.0,
    ('T', 'W'): 2.2,  ('W', 'T'): 2.2,
    ('H', 'W'): -0.6, ('W', 'H'): -0.6,
    ('H', 'T'): 0.0,  ('T', 'H'): 0.0,
}


def pair_energy(a, b):
    return EPS[(a, b)]


class System:
    def __init__(self, L, n_chains, n_head, n_tail, rng=None):
        self.L = L
        self.rng = rng or random.Random(12345)
        self.chain_len = n_head + n_tail
        self.type_seq = ['H'] * n_head + ['T'] * n_tail  # head block then tail block
        self.grid = np.full((L, L), 'W', dtype='<U1')
        self.chains = []  # list of lists of (x,y), head-to-tail order
        self._place_chains(n_chains)

    def _neighbors(self, x, y):
        L = self.L
        return ((x + 1) % L, y), ((x - 1) % L, y), (x, (y + 1) % L), (x, (y - 1) % L)

    def _interaction_neighbors(self, x, y):
        """Orthogonal neighbors at full weight + diagonal neighbors at half
        weight. Reduces the square-lattice bias toward axis-aligned, faceted
        aggregate shapes by giving beads a more isotropic local environment
        (closer to a real continuum coordination shell) without changing
        chain bond connectivity itself."""
        L = self.L
        ortho = [((x + 1) % L, y, 1.0), ((x - 1) % L, y, 1.0),
                 (x, (y + 1) % L, 1.0), (x, (y - 1) % L, 1.0)]
        diag = [((x + 1) % L, (y + 1) % L, 0.5), ((x + 1) % L, (y - 1) % L, 0.5),
                ((x - 1) % L, (y + 1) % L, 0.5), ((x - 1) % L, (y - 1) % L, 0.5)]
        return ortho + diag

    def _place_chains(self, n_chains):
        placed = 0
        attempts = 0
        while placed < n_chains and attempts < n_chains * 500:
            attempts += 1
            x0, y0 = self.rng.randrange(self.L), self.rng.randrange(self.L)
            if self.grid[x0, y0] != 'W':
                continue
            path = [(x0, y0)]
            visited = {(x0, y0)}
            ok = True
            for _ in range(self.chain_len - 1):
                x, y = path[-1]
                opts = [p for p in self._neighbors(x, y)
                        if p not in visited and self.grid[p] == 'W']
                if not opts:
                    ok = False
                    break
                nxt = self.rng.choice(opts)
                path.append(nxt)
                visited.add(nxt)
            if not ok:
                continue
            for (x, y), t in zip(path, self.type_seq):
                self.grid[x, y] = t
            self.chains.append(path)
            placed += 1
        if placed < n_chains:
            raise RuntimeError(f"Only placed {placed}/{n_chains} chains -- lattice too crowded")

    def _delta_E(self, changes):
        # changes: dict site -> (old_type, new_type)
        edges = set()  # (site_a, site_b, weight)
        for (x, y) in changes:
            for (nx, ny, w) in self._interaction_neighbors(x, y):
                a, b = (x, y), (nx, ny)
                key = (a, b, w) if a < b else (b, a, w)
                edges.add(key)
        dE = 0.0
        for (p, q, w) in edges:
            old_p = changes[p][0] if p in changes else self.grid[p]
            new_p = changes[p][1] if p in changes else self.grid[p]
            old_q = changes[q][0] if q in changes else self.grid[q]
            new_q = changes[q][1] if q in changes else self.grid[q]
            dE += w * (pair_energy(new_p, new_q) - pair_energy(old_p, old_q))
        return dE

    def _apply(self, changes):
        for (x, y), (_, new_t) in changes.items():
            self.grid[x, y] = new_t

    def try_end_move(self, ci):
        chain = self.chains[ci]
        end_is_head = self.rng.random() < 0.5
        idx = 0 if end_is_head else len(chain) - 1
        anchor_idx = 1 if end_is_head else len(chain) - 2
        end_site = chain[idx]
        anchor_site = chain[anchor_idx]
        opts = [p for p in self._neighbors(*anchor_site)
                if p != end_site and self.grid[p] == 'W']
        if not opts:
            return
        new_site = self.rng.choice(opts)
        t = self.type_seq[idx]
        changes = {end_site: (t, 'W'), new_site: ('W', t)}
        dE = self._delta_E(changes)
        if dE <= 0 or self.rng.random() < np.exp(-dE):
            self._apply(changes)
            chain[idx] = new_site

    def try_kink_move(self, ci):
        chain = self.chains[ci]
        n = len(chain)
        if n < 3:
            return
        i = self.rng.randrange(1, n - 1)
        prev_site, cur_site, next_site = chain[i - 1], chain[i], chain[i + 1]
        px, py = prev_site
        nx_, ny_ = next_site
        # only a kink if prev and next are diagonal (not colinear, not same point)
        if px == nx_ or py == ny_:
            return
        cand1 = (px, ny_)
        cand2 = (nx_, py)
        cand = cand1 if cand1 != cur_site else cand2
        if self.grid[cand] != 'W':
            return
        t = self.type_seq[i]
        changes = {cur_site: (t, 'W'), cand: ('W', t)}
        dE = self._delta_E(changes)
        if dE <= 0 or self.rng.random() < np.exp(-dE):
            self._apply(changes)
            chain[i] = cand

    def try_rigid_translation(self, ci):
        chain = self.chains[ci]
        own_sites = set(chain)
        d = self.rng.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
        L = self.L
        new_positions = [((x + d[0]) % L, (y + d[1]) % L) for (x, y) in chain]
        for p in new_positions:
            if p not in own_sites and self.grid[p] != 'W':
                return
        old_types = {s: t for s, t in zip(chain, self.type_seq)}
        new_types = {s: t for s, t in zip(new_positions, self.type_seq)}
        all_sites = set(chain) | set(new_positions)
        changes = {}
        for s in all_sites:
            old_t = old_types[s] if s in old_types else self.grid[s]
            new_t = new_types[s] if s in new_types else 'W'
            changes[s] = (old_t, new_t)
        dE = self._delta_E(changes)
        if dE <= 0 or self.rng.random() < np.exp(-dE):
            self._apply(changes)
            self.chains[ci] = new_positions

    def sweep(self):
        n = len(self.chains)
        order = list(range(n))
        self.rng.shuffle(order)
        for ci in order:
            self.try_end_move(ci)
            self.try_kink_move(ci)
            self.try_kink_move(ci)
            if self.rng.random() < 0.15:
                self.try_rigid_translation(ci)

    def run(self, n_sweeps):
        for _ in range(n_sweeps):
            self.sweep()

    def total_energy(self):
        E = 0.0
        L = self.L
        for x in range(L):
            for y in range(L):
                t1 = self.grid[x, y]
                # each orthogonal edge once (right, up), each diagonal edge once (down-right, up-right)
                for (nx, ny, w) in [((x + 1) % L, y, 1.0), (x, (y + 1) % L, 1.0),
                                     ((x + 1) % L, (y + 1) % L, 0.5), ((x + 1) % L, (y - 1) % L, 0.5)]:
                    E += w * pair_energy(t1, self.grid[nx, ny])
        return E
