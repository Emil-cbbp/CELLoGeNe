# CELLoGeNe
The CELLoGeNe (Computation of Energy Landscapes of Logical Gene Networks) software, written and developed by Emil Andersson* and Mattias Sjö**

\*Computational Biology and Biological Physics group, Department of Astronomy and Theoretical Physics, Lund University.
** Theoretical Particle Physics group, Department of Astronomy and Theoretical Physics, Lund University.

DOI: https://doi.org/10.1016/j.isci.2022.104743

CELLoGeNe consists of three parts:
  * Calculating energy landscapes from a given gene regulatory network (GRN). CELLoGeNe can either calculate the energy landscape for a single specified configuration of logical operators, or test different configurations either exhaustively or a set of random configurations. (Simulation.py, Network.py, RunData.py, Operators.py, MatrixGenerator.py, lgnReader.py, GeneTree.py, ExpressionSpace.py)
  * Constructing multidimensional visualisation plots of discrete energy landscapes. (GraphPlotter.py)
  * Perform simmulations of cells moving in the energy landscape under influence of noise. (CellSimulation.py, Simulate_cells_demo.py)

Demos for each of the three parts are included. The GRN and other details are specified in a .lgn file. Here a demo of a toy GRN is included in HowTo.lgn. Run the demo from the command line with:

python Simulation.py HowTo.lgn

Create a plot for an energy landscape with:

python Graphplotter.py

The resulting file can be compiled with pdfLaTeX

Run the cell simulation demo with:

python Simulate_cells_demo.py
