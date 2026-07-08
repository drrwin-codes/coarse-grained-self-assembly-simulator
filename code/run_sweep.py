import sys, time, csv, random
sys.path.insert(0, '/home/claude/self_assembly')
from lattice_mc import System
from cluster_analysis import classify_system_with_vesicles
from collections import Counter
import statistics

L = 40
CHAIN_LEN = 8
SWEEPS = 6000
SEEDS = [11, 22, 33, 44]

g_configs = [(7, 1), (6, 2), (5, 3), (4, 4), (3, 5), (2, 6), (1, 7)]  # (NH, NT)

n_chains = int(sys.argv[1])  # pass one concentration level per invocation
phi = n_chains * CHAIN_LEN / (L * L)
outpath = f'/home/claude/self_assembly/sweep_chunk_{n_chains}.csv'

fieldnames = ['phi', 'g', 'nh', 'nt', 'mean_agg_number', 'std_agg_number',
              'mean_max_agg_number', 'majority_label', 'label_agreement',
              'vesicle_fraction', 'per_replicate_dominant']

t_start = time.time()
with open(outpath, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for nh, nt in g_configs:
        g = nt / (nh + nt)
        rep_mean_agg, rep_max_agg, rep_dominant, rep_has_vesicle = [], [], [], []
        for seed in SEEDS:
            s = System(L=L, n_chains=n_chains, n_head=nh, n_tail=nt,
                       rng=random.Random(seed))
            s.run(SWEEPS)
            infos, vesicle_info = classify_system_with_vesicles(s)
            total_chains = sum(r['agg_number'] for r in infos)
            chains_by_label = {}
            for r in infos:
                chains_by_label[r['label']] = chains_by_label.get(r['label'], 0) + r['agg_number']
            dominant_label = max(chains_by_label, key=chains_by_label.get)
            mean_agg_number = sum(r['agg_number'] ** 2 for r in infos) / total_chains
            max_agg_number = max(r['agg_number'] for r in infos)

            rep_mean_agg.append(mean_agg_number)
            rep_max_agg.append(max_agg_number)
            rep_dominant.append(dominant_label)
            rep_has_vesicle.append(len(vesicle_info) > 0)

        mean_of_means = statistics.mean(rep_mean_agg)
        std_of_means = statistics.stdev(rep_mean_agg) if len(rep_mean_agg) > 1 else 0.0
        mean_of_max = statistics.mean(rep_max_agg)
        vesicle_fraction = sum(rep_has_vesicle) / len(rep_has_vesicle)
        majority_label = Counter(rep_dominant).most_common(1)[0][0]
        label_agreement = Counter(rep_dominant).most_common(1)[0][1] / len(rep_dominant)

        row = {
            'phi': phi, 'g': g, 'nh': nh, 'nt': nt,
            'mean_agg_number': mean_of_means, 'std_agg_number': std_of_means,
            'mean_max_agg_number': mean_of_max,
            'majority_label': majority_label, 'label_agreement': label_agreement,
            'vesicle_fraction': vesicle_fraction,
            'per_replicate_dominant': ';'.join(rep_dominant),
        }
        writer.writerow(row)
        f.flush()
        print(f"phi={phi:.3f} g={g:.3f}  mean_agg#={mean_of_means:.1f}+/-{std_of_means:.1f}  "
              f"mean_max#={mean_of_max:.1f}  majority={majority_label} "
              f"(agree={label_agreement:.2f})  vesicle_frac={vesicle_fraction:.2f}  "
              f"elapsed={time.time()-t_start:.0f}s")

print(f"\nChunk phi={phi:.3f} done in {time.time()-t_start:.0f}s. Saved {outpath}")
