# Winner-take-all networks for graph coloring (wta-graphcoloring)

## Overview
This directory contains Matlab code to simulate solving constraint satisfaction problems with networks of winner-take-all networks. This work is the result of a collaboration between [Ueli Rutishauser](http://rutishauserlab.org), [Jean-Jacques Slotine](http://web.mit.edu/nsl/www/), and [Rodney Douglas](http://www.ini.uzh.ch).

The purpose of this released code is to reproduces key aspects of this paper:
>Solving constraint-satisfaction problems with distributed neocortical-like neuronal networks. Rutishauser U, Jean-Jacques Slotine, Rodney Douglas Neural computation, 30(5):1359-1393 (2018). [Fulltext](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5930080/)

The key theoretical concepts used to make these networks work are explained in this paper:
>Computation in dynamically bounded asymmetric systems. U. Rutishauser, JJ. Slotine, RJ. Douglas. PLOS Computational Biology, 11(1): e1004039 (2015). [Fulltext](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004039)

Also see this page for [code related to our earlier WTA papers](https://services.ini.uzh.ch/~urut/DFAWTA/).

This code is released as-is for academic purposes. Feel free to use and re-use as you see fit.
Uploaded January 2019, Ueli Rutishauser. www.rutishauserlab.org 

## How to get started
Download the entire repository. Then, adjust and run setpath_graph.m to set the path. Then, execute one of the main functions listed below.

## List of key functions (main)

```mainSimpleRecurrence.m``` Simple example that illustrates how two coupled WTAs can hold a state after offset of input.

```simpleRecurrence_DBcells_twoWinnersOnly_Fig2.m``` Simple example that illustrates the concept of negative constraint cells (also called "DB cells" in the code).
This function reproduces Figure 2 of the paper.

```mainWTA_graphColoring_Fig3.m``` Simple graph coloring example (4 nodes). Reproduces Figure 3 of the paper.

```mainWTA_graphColoring_Fig4.m``` 4-coloring of random planar graphs of arbitrary size. Reproduces Figure 4 of the paper.

```mainWTA_graphColoring_Fig5.m``` Solving the MIS problem with graph coloring. Reproduces Figure 5 of the paper.

```mainWTA_graphColoring_Fig6.m``` Solves sudoku problems. Reproduces Figure 6 of the paper.

## List of key functions (internal, called from others)

```mainWTA_graphColoring_run1.m``` The main function. Setups up the network, executes it (using runN_*, see below), and evaluates the result.

```runN_WTA_withDBs_optimized_shunt.m``` Function that does the numerical integration of the network.

## 3rd Party tools used
Apart from Matlab, the following 3rd Party tools are used.

Required (included in release for convenience):
1. [Boost Graph Library, MatlabBGL version] (http://dgleich.github.io/matlab-bgl/)
2. [GraphViz2Mat](https://www.mathworks.com/matlabcentral/fileexchange/4518-matlab-graphviz-interface)

Optimal:
1. JFLAP (+Java) (http://jflag.org)
2. graphviz

## XML Graph format

The WTA code stores and reads arbitrarily generated random planar graphs (or the sudoku graph) stores as JFF files.
JFF files are XML files that describe a graph. A convenient way to edit and view these graphs is JFLAP, a java utility.

Graph files that are used: 
1. sudoku_s1.jff
2. maxIndpSet_Fig5A.jff  (Fig 5A graph)
3. graph3.jff (Fig 3A graph)
4. exportRand_*.jff  (Fig 4 and 5 randomly generated planar graphs)

## Layouting of graphs

Two ways to layout (plot) the resulting graph coloring solutions are implemented:
Using the boost graph library (matlab_bgl) or by exporting a .dot file and running graphviz (dot).

