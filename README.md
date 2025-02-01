# Planar Maps Library for SageMath

This project is a **SageMath library** for manipulating **planar maps**.  Maps are **embeddings of graphs on surfaces**, commonly used in **graph theory, combinatorics, and topology**.  
This library provides **tools for creating, visualizing, and manipulating planar maps** within the SageMath environment.

---

## Documentation

The implementation of planar maps is primarily based on the following references:

- **[1]** *Handbook on Enumerative Combinatorics*, Gilles Schaeffer, 2015 [Website](https://www.taylorfrancis.com/chapters/edit/10.1201/b18255-11/planar-maps-gillesschaeffer)
- **[2]** *An Introduction to Map Enumeration*, Guillaume Chapuy, 2014 [Website](https://www.irif.fr/~chapuy/AEC/coursHagenberg.pdf)
- **[3]** *Enumeration and Random Generation of Planar Maps*, Ivan Geffner, Master’s Final Thesis, Universitat Politècnica de Catalunya, 2014 [Website](https://upcommons.upc.edu/handle/2099.1/23114)

---

## Getting Started

This library requires a **working SageMath installation**.  
Before using it, we highly recommend reading the official SageMath documentation:

- **[SageMath README](https://doc.sagemath.org/html/en/README.html)**
- **[Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html)**

---

## Overview

The library contains **four main classes**:

- [**Labelled Map**](#labelled-map)
- [**Rooted Map**](#rooted-map)
- [**Mutable Labelled Map**](#mutable-labelled-map)
- [**Map Generator**](#map-generator)

Each class provides specific functionalities for working with planar maps.

---

## Labelled Map

A **Labelled Map** is a **graph representation** equipped with a **rotation system** and an arbitrary **labelling of its \(2n\) half-edges** with numbers in `[1 ... 2n]`.  

This class provides **two constructors**:

1. **From permutations** (`σ` and `α`), following the notation in **[2]**.
2. **From an adjacency list**.

### **Main Methods**
- Compute:
  - **Face count**, **number of nodes**, **edges**, **genus**, **diameter**.
- Construct various **derived maps**:
  - **Spanning Tree**, **Dual Graph**, **Derived Map**, **Incidence Map**.
- **Visualization**:
  - `show()` method for **plotting the map**.

---

## Rooted Map

A **Rooted Map** is an **equivalence class of labelled maps** under relabelling of `[1 ... 2n]`, preserving the **root half-edge** (labelled as `1`).  
The **half-edge labelled `1`** serves as the **root** of the map.

This class **inherits** from `LabelledMap`.

---

## Mutable Labelled Map

This class **inherits** from `LabelledMap` and provides **methods for modifying** the map, such as:

- **Adding or removing edges**.
- **Deleting vertices**.
- **Contracting edges**.
- **Contracting faces**.

---

## Map Generator

Provides methods for **generating maps**, including:

- A method for generating a **uniformly random rooted map** with a fixed number of edges.

---

## Contact

For **questions, issues, or feature requests**, please **open an issue** on the **GitHub repository**.
