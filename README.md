
[![Delaunay triangulation](images/lifted.png?raw=true "Delaunay triangulation")]
===============================================================================

This program shows a visualization of computation of Delaunay triangulation
by computing the convex hull of lifted 2D points to a 3D paraboloid. 

Qhull http://www.qhull.org/ is used to do most of the heavy lifting stuff. 

Part of the code is from Computational Geometry in C by 

There is no input argument. You can modify the visualization via the following
variables in main.c 

const int width  = 1;
const int height = 1;
const int sample = 1000;

Here are the keys that you can use to control the visualization

GUI keys:
r:   select a random circumcircle
a:   show/hide all circumcircles
v:   show Voronoi verices
?:   show this message
esc: quit

===============================================================================

Let me know if there are problems/bug.
jmlien@cs.gmu.edu
