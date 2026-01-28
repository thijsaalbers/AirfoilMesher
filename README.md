# Airfoil Mesher

A simple python script capable of automatically creating meshes of basic airfoil geometries.

Note how more advanced, and more efficient codes with the same purpose already exist. This code is purely a means for myself to practice my coding skills as well as create a really basic utility for my own purposes.

---

## Required libraries

The required python libraries are as follows

```bash
    pip install gmsh
    pip install matplotlib
    pip install numpy
```

---

## User guide

The codes main functionality is creating meshes of airfoil geometries. For now this is limited to the NACA 4 digit series, more airfoil geometries might be added later if required. 

The meshing of the airfoil geometry currently has two options

1. Creating unstructured grids using GMSH automatic mesh generation
2. Creating an O-grid using elliptical grid generation

To call the code and create an unstructured grid of a NACA-2412 airfoil the following command can be used

```bash
    python3 AirfoilMesher.py NACA2412 gmsh
```

To create a structured O-grid of the same geometry, the following command can be used

```bash
    python3 AirfoilMesher.py NACA2412 elliptic
```

To access specific meshing options such as the element size or the meshing algorithm, please see the corresponding python code files themselves. For future work this might be added as command line arguments.

For both meshing options, the GMSH FLTK GUI is opened upon finishing mesh generation. In this way, additional mesh manipulation through this GUI is available as well as saving the completed mesh to different output formats courtesy of GMSH.

---



