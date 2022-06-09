# ap_data_generation

This repository contains the data generation pipeline for the preprint [Phase Object Reconstruction for 4D-STEM using Deep Learning](https://arxiv.org/abs/2202.12611).
It depends on the repositories [multem](https://github.com/Ivanlh20/multem) and [atomic_specimen_creation](https://github.com/ThFriedrich/atomic_specimen_creation)
The atomic structures are created from cif files, which were generated with [pymatgen](https://pymatgen.org/). The cif files used are provided in this repository too via git lfs.

To clone with submodules run:
```bash
git clone --recurse-submodules -j3 git@github.com:ThFriedrich/ap_data_generation.git
```
