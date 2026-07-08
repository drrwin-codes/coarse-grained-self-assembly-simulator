# Coarse Grained Self-Assembly Simulation

**r_q: can a lightweight, from-scratch simulation reproduce the known phase map of these structures, using nothing but molecular geometry as an input?**

The spontaneous organization of **`amphiphilic molecules`** into well-defined macroscopic topologies, such as **`spherical micelles`, `wormlike cylinders`**, and **`closed bilayer vesicles`**, outlines essential biological motifs and industrial formulations. 

> While high-fidelity molecular dynamics simulations offer exceptional atomistic detail, their extensive computational demands mask the primary physical invariants governing these phase transitions. Here, we present an ultra-lightweight, `2D lattice Monte Carlo simulator` designed to establish whether complex, multi-state phase landscapes can be accurately recovered solely from local molecular geometry and concentration.

By representing amphiphiles as minimal bead-chains with simple solvophobic energy rules, our model successfully captures the emergence of diverse self-assembled aggregates without imposing explicit morphological templates. We demonstrate that this minimal engine quantitatively matches the analytical boundaries predicted by `Israelachvili’s classical packing parameter framework`.

> To facilitate broader scientific exploration and educational engagement, we integrate this discrete physical engine into an optimized, browser-based real-time visualization platform. Our findings prove that local steric constraints and interfacial energy minimization are completely sufficient to dictate global amphiphilic phase topographies, isolating the essential noise from the structural signal in self-assembly phenomena


