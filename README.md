# HeatDiffusion
This project is based on an ACM Transactions on Graphics article by Crane et al., [Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow](https://dl.acm.org/citation.cfm?doid=2516971.2516977) . Keenan Crane introduced the heat method for computing the deodesic distance efficiently. In practice, this method updates the distance an order of magnitude faster than the previously existing methods. This project is a reproduction of the algorithm proposed in the paper.
 <img src="/demo/heatdemo.gif">

## Instruction
The project uses Java and the JOGL openhl bindings.
Download the code and add the jogl, vecmath, and mintools jars to the root level of the project folder.
The utility files came from the 2.3.2 release folder (https://jogamp.org/wiki/index.php/Release_2.3.2) on the jogamp site.  The two relevant files are 7unzipped and attached [here](https://github.com/dorhelium/HeatDiffusion/tree/master/utilities), so you don't need to install 7zip if you don't have it).  The fat jar contains all the necessary dlls and shared objects to work on all common platforms.

## Contents

### MeshProcessingApp
Contains the main function and creates the view.
### HalfEdge 
Half edge class including code to draw the half edge to help with debugging.
### Face 
Face class with additional members to store the center position, the area, and the gradient of the heat value across the mesh.
### Vertex
Vertex class with additional members for storing a "smooth" surface normal, Laplacian coefficients, heat values, initial conditions, etc.
### PolygonSoup 
Simple parser for obj files, scale and center vertices after loading.
### HEDS 
The half edge data structure (to construct from a PolygonSoup), which also contains a number of method stubs for the geometry processing.
### MeshDrawHeatGeo 
This class sets up native IO vertex and index buffers for drawing, loads and links the pflHeatGeo vertex and fragment program (fpl means per fragment lighting), and sets up a collection of uniform values (e.g., model view and projection matrices) and attribute buffers for drawing.
### MeshDrawPicking 
This class sets up a picking GLSL program, and likewise sets up buffers to draw every triangle in a unique colour to allow selection of mesh faces. While vertex picking is needed, the barycentric coordinates of a triangle allows for the closest vertex to the click point to be selected.


    
