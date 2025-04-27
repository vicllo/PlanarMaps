# Maps Library for SageMath

This project is a SageMath library for manipulating maps. Maps are embeddings of graphs on surfaces, commonly used in graph theory, combinatorics, and topology.  
This library provides tools for creating, visualizing, and manipulating planar maps within the SageMath environment.

---

## Documentation

The implementation of maps is primarily based on the following references:

- [1] *Handbook on Enumerative Combinatorics*, Gilles Schaeffer, 2015 [Website](https://www.taylorfrancis.com/chapters/edit/10.1201/b18255-11/planar-maps-gillesschaeffer)

- [2] *An Introduction to Map Enumeration*, Guillaume Chapuy, 2014 [Website](https://www.irif.fr/~chapuy/AEC/coursHagenberg.pdf)

- [3] *Enumeration and Random Generation of Planar Maps*, Ivan Geffner, Master’s Final Thesis, Universitat Politècnica de Catalunya, 2014 [Website](https://upcommons.upc.edu/handle/2099.1/23114)

- [4] *Optimal Coding and Sampling of Triangulations*, Dominique Poulalhon and Gilles Schaeffer, LIX – CNRS Ecole Polytechnique [Website](https://www.lix.polytechnique.fr/Labo/Dominique.Poulalhon/Articles/PoSc_cod_ICALP.pdf)

---

## Getting Started

This library requires a working SageMath installation.  
Before using it, we highly recommend reading the official SageMath documentation:

- [SageMath README](https://doc.sagemath.org/html/en/README.html)
- [Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html)

We recommend using Conda Forge:

- [Sage Installation Guide (Conda Forge)](https://doc.sagemath.org/html/en/installation/conda.html)

After completing the installation, check the file `example.py` for usage examples or run it by : python example.py to see them in action.

---

## Overview

The library contains 6 main classes:

- [LabelledMap](#labelled-map)
- [RootedMap](#rooted-map)
- [MutableLabelledMap](#mutable-labelled-map)
- [PrimitiveMutableLabelledMap](#primitive-mutable-labelled-map)
- [MapGenerator](#map-generator)
- [DynamicShow](#dynamicshow)

Each class provides specific functionalities for working with maps.

---

## Labelled Map

A *Labelled Map* is a graph representation equipped with a rotation system and an arbitrary labeling of its \(2n\) half-edges(also called demi-edges) with numbers in `[1 ... 2n]`.

This class provides two constructors:

1. From permutations (`σ` and `α`), following the notation in [2].
2. From an adjacency list.

### Main Methods

- Compute:
  - Face count, number of nodes, edges, genus, diameter, etc.
- Construct various derived maps:
  - Spanning Tree, Dual Graph, Derived Map, Incidence Map, etc.
- Visualization:
  - `show()` method for plotting the map.

It also provides a more user-friendly way to interact with demi-edges other than raw indices, by allowing the use of `TopologicalDemiEdge`.

---

## Rooted Map

A *Rooted Map* is an equivalence class of labelled maps under relabeling of `[1 ... 2n]`, preserving the root half-edge (labelled as `1`).  
The half-edge labelled `1` serves as the root of the map. All methods that return a map will return a `RootedMap`.

This class inherits from `LabelledMap`.

---

## Mutable Labelled Map

This class inherits from `LabelledMap` and provides methods for modifying the map, such as:

- Adding or removing edges
- Deleting vertices
- Contracting edges
- Contracting faces
- etc.

Some notes:

- All methods that return a map will return a `LabelledMap`, not a `MutableLabelledMap`.
- Most methods that modify the map will alter its labels. To keep track of demi-edges even when labels change, please use `MutableTopologicalDemiEdge`, which guarantees that the label associated with it is updated accordingly.

---


## Primitive Mutable Labelled Map

This class inherits from `LabelledMap` and provides methods for modifying the map, it is much more primitive than MutableLabelledMap and place more responsability on the user and if you aren't careful the map may become not stable , some methods of LabelledMap are not implemented for this class and will raise an error.Only useful compared to MutableLabelledMap because on some operations it is O(1) instead of O(log(m)), if it isn't critical for your work we advised you to use MutableLabelledMap instead.

Some notes:

- All methods that return a map will return a `LabelledMap`, not a `PrimitiveMutableLabelledMap`.
- Most methods that modify the map will alter its labels. To keep track of demi-edges even when labels change, please use `PrimitiveMutableTopologicalDemiEdge`, which guarantees that the label associated with it is updated accordingly.



---

## Map Generator

Provides methods for generating maps or map-related objects, including:

- A method for generating  uniformly random rooted map with a fixed number of edges.
- A method for generating  uniformly random rooted tree of fixed size.
- A method for generating  uniformly random Triangulations of fixed size  
- etc.

---

## DynamicShow

This class provides a more advanced and customizable way for users to display maps than the default `show` method.

---

## Other Classes

- `TopologicalDemiEdge`: A more user-friendly way to interact with demi-edges than using raw indices.
- `MutableTopologicalDemiEdge`: Inherits from `TopologicalDemiEdge`; adds methods only possible with mutable maps.
- `MapPermutation`: A custom implementation of permutations used internally by the library.
- `RotatingPermutation`: Inherits from `MapPermutation`; a special kind of mutable permutation with performant query time \(O(\log n)\) operations (e.g., checking if two indices are on the same cycle) and efficient operations like deleting or adding indices in a cycle.

- `PrimitiveRotatingPermutation`: Inherits from `MapPermutation`; a more primitive version of RotatingPermutation it support modification operations in O(1).

- `CustomSwap`: Inherits from `MapPermutation`; a special class for transpositions, more efficient in some cases than using `MapPermutation` directly.
- `PermutationUtilsAbstractor`: An internal class used to abstract some methods related to map permutations.
- `RotatingPermutationUtilsAbstractor`: Inherits from `PermutationUtilsAbstractor`, adapted for `RotatingPermutation`.
- `PrimitiveRotatingPermutationUtilsAbstractor`: Inherits from `PermutationUtilsAbstractor`, adapted for `PrimitiveRotatingPermutation`.

- `CyclicChainedList`: An implementation of a mutable cyclic list used in `RotatingPermutation`.
- `CyclicUtilsProvider`: A class that abstracts operations (e.g., querying whether two indices are on the same cycle), used in `RotatingPermutation`.
- `SplayTree`: An implementation of a modified [splay tree](https://en.wikipedia.org/wiki/Splay_tree) data structure with some tweaks useful in `CyclicUtilsProvider`.

---

## Contact

For questions, issues, or feature requests, please open an issue on the GitHub repository.

