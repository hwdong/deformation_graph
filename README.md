# deformation_graph
This is my deformation graph code for the implementation of Embedded deformation. The deformation graph will look like the following picture:

![alt tag](https://github.com/hwdong/deformation_graph/blob/master/geo_sample.jpg) 

The deformation graph is generated using farthest point sampling algorithm and 
the code is based on the fast marching code from Professor Gabriel Peyre to compute geodesic distance of a triangular mesh.  

Please look at the gen_deformation_graph.h on how to generate defoormation graph. You need not use the test_FPS.cpp which is just for test fast marching algorithm (farthest point sampling algorithm).


Reference

1. Embedded Deformation for Shape Manipulation. 
    https://graphics.ethz.ch/~sumnerb/research/embdef/Sumner2007EDF.pdf
    
2. Geodesic Computations on 3D Meshes.http://www.cmap.polytechnique.fr/~peyre/geodesic_computations/