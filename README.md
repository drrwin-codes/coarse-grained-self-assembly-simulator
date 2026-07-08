# Coarse Grained Self-Assembly Simulation

**r_q: can a lightweight, from-scratch simulation reproduce the known phase map of these structures, using nothing but molecular geometry as an input?**

`Amphiphilic molecules` don't need to be told how to organize. Left in solution, they spontaneously sort themselves into spherical `micelles`, wormlike `cylinders`, or closed `bilayer vesicles`, structures that underpin everything from cell membranes to laundry detergent.

`Molecular dynamics simulations` can capture this `self-assembly` in exquisite atomistic detail, but that detail comes at a cost: 

> **the computational overhead often buries the simple physical rules actually driving the transitions between these phases.**

**I wanted to strip the problem down and ask a more basic question:**

> **how much of this behavior can be recovered from geometry and concentration alone, with none of the chemical bookkeeping?**

To find out, I built a `minimal 2D lattice Monte Carlo simulator`. Amphiphiles are represented as short `bead-chains` governed by simple `solvophobic energy rules`, nothing more. No morphology is hardcoded into the model. And yet the full range of aggregate shapes emerges on its own, matching the phase boundaries predicted by `Israelachvili's classical packing parameter framework` with surprising fidelity. I also built a lightweight, browser-based visualization layer on top of the engine, partly for my own exploration, partly so others can watch these phases form in real time without needing a simulation background.

---

> <img width="2880" height="1627" alt="image" src="https://github.com/user-attachments/assets/43b811d7-8272-4b3c-885b-06ecc5768e8a" />

---

 `The upshot`: **local packing constraints and interfacial energy minimization, on their own, are enough to determine the global phase landscape of amphiphile self-assembly. 

---
> Phase behavior isn't an emergent property that requires the full weight of atomistic simulation to predict; a packing parameter and a concentration are, in this minimal picture, doing most of the explanatory work; the chemical identity of the tail, the exact solvent composition, the fine details MD captures so well, seems to matter less for determining which phase you land in, and more for the finer details of *how* you get there.
 ---
There's something clarifying about watching a system this stripped-down still find its way to spheres, cylinders, and bilayers. Hope you find it cool too :)

