# HeatDiffusion

### Instruction
Download the code and add the jogl, vecmath, and mintools jars to the project.

### Contents
MeshProcessingApp - Contains the main function and creates the view.
HalfEdge - Half edge class, much like what was seen in the lectures, and including code to draw the half edge to help with debugging. Feel free to add members as you see fit (e.g., cache important values for expensive computations).
Face - Face class, much like what we have seen in the lectures (e.g., having an example half edge that has this face as its leftface), but with additional members to store the center position, the area, and the gradient of the heat value across the mesh.
Vertex - Vertex class, much like what we have seen in the lectures (e.g., having an example half edge that has this vertex as its head), but with additional members for storing a "smooth" surface normal, Laplacian coefficients, heat values, initial conditions, etc.
PolygonSoup - Simple parser for obj files, which you will modify to scale and center vertices after loading.
HEDS - The half edge data structure (to construct from a PolygonSoup), which also contains a number of method stubs for you to complete the geometry processing portion of the assignment objectives..
MeshDrawHeatGeo - This class sets up native IO vertex and index buffers for drawing, loads and links the pflHeatGeo vertex and fragment program (fpl here means per fragment lighting, which is an objective below), and sets up a collection of uniform values (e.g., model view and projection matrices) and attribute buffers for drawing.
MeshDrawPicking - Much like the MeshDrawHeatGeo class, this class sets up a picking GLSL program, and likewise sets up buffers to draw every triangle in a unique colour to allow selection of mesh faces. While vertex picking is needed in the assignment, the barycentric coordinates of a triangle allows for the closest vertex to the click point to be selected.
