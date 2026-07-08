import sys, time
sys.path.insert(0, '/home/claude/self_assembly')
from lattice_mc import System
import numpy as np
import matplotlib.pyplot as plt

L = 40
N_CHAINS = 45
CHAIN_LEN = 8
SWEEPS = 6000

cases = [
    ("head_heavy_g025", 6, 2),  # NH=6, NT=2 -> g = NT/(NH+NT) = 0.25
    ("tail_heavy_g075", 2, 6),  # NH=2, NT=6 -> g = 0.75
]

fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))

for ax, (label, nh, nt) in zip(axes, cases):
    g = nt / (nh + nt)
    print(f"Running {label} (g={g:.2f}) ...")
    sysm = System(L=L, n_chains=N_CHAINS, n_head=nh, n_tail=nt, )
    e0 = sysm.total_energy()
    t0 = time.time()
    sysm.run(SWEEPS)
    e1 = sysm.total_energy()
    dt = time.time() - t0
    print(f"  E: {e0:.1f} -> {e1:.1f}  ({SWEEPS} sweeps, {dt:.1f}s)")

    img = np.zeros((L, L, 3))
    color = {'W': (1, 1, 1), 'H': (0.25, 0.45, 0.95), 'T': (0.85, 0.2, 0.2)}
    for x in range(L):
        for y in range(L):
            img[y, x] = color[sysm.grid[x, y]]  # transpose for imshow orientation

    ax.imshow(img, origin='lower', interpolation='nearest')
    ax.set_title(f"g = {g:.2f}  (NH={nh}, NT={nt})\nfinal energy = {e1:.0f}")
    ax.set_xticks([])
    ax.set_yticks([])

axes[0].set_ylabel("head-heavy geometry")
fig.suptitle(f"2D lattice MC self-assembly after {SWEEPS} sweeps  "
             f"(blue=head, red=tail, white=solvent)")
plt.tight_layout()
plt.savefig('/home/claude/self_assembly/demo_snapshots.png', dpi=150)
print("saved demo_snapshots.png")
