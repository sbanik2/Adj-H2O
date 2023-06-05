# AdjH~2~O

## Description


Python code for calculating the adjacency matrix of H2O, with and without considering the hydrogen bond in the atomic configuration, using the criteria defined in the [Paper](https://doi.org/10.1038/379055a0).

The code takes the atomic configurations as a [Pymatgen](https://pymatgen.org/)  structure


 ``` python
 
from adjacency import  get_adjacency
from pymatgen.core import Structure



structure = Structure.from_file("POSCAR") 

A = get_adjacency(
                 structure,
                 rcut_OH=1.1,
                 rcut_OO=3.5,
                 phi = 30,
                 consider_Hbond=True
                 
                    )
 
 
 ```



## Reference

Luzar, A., & Chandler, D. (1996). Hydrogen-bond kinetics in liquid water. Nature, 379(6560), 55-57.
