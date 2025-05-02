// Geometry parameters
L = 0.8;       // Length of domain
H = 0.41;      // Height of domain
r = 0.05;      // Cylinder radius
xc = 0.2;      // Cylinder center x
yc = 0.2;      // Cylinder center y

// Mesh resolution
res = 0.01;

// Outer rectangle
Point(1) = {0, 0, 0, res};
Point(2) = {L, 0, 0, res};
Point(3) = {L, H, 0, res};
Point(4) = {0, H, 0, res};

// Cylinder center (used as arc center)
Point(5) = {xc, yc, 0, res};

// Points on the circle (4 quadrants)
Point(6) = {xc + r, yc, 0, res};
Point(7) = {xc, yc + r, 0, res};
Point(8) = {xc - r, yc, 0, res};
Point(9) = {xc, yc - r, 0, res};

// Arcs forming the circle
Circle(10) = {6, 5, 7};
Circle(11) = {7, 5, 8};
Circle(12) = {8, 5, 9};
Circle(13) = {9, 5, 6};

// Outer boundary lines
Line(14) = {1, 2};
Line(15) = {2, 3};
Line(16) = {3, 4};
Line(17) = {4, 1};

// Create line loops
Line Loop(18) = {14, 15, 16, 17};
Line Loop(19) = {10, 11, 12, 13};

// Define the surface (domain minus cylinder)
Plane Surface(20) = {18, 19};

// Recombine to quadrilaterals
Recombine Surface{20};

// Optional: create mild non-orthogonality via displacement field
// Field[1] = MathEval;
// Field[1].F = "0.01*sin(8*x)*cos(5*y)";
// DisplacementField = 1;

// Physical groups (optional for boundary conditions)
Physical Line("inlet") = {17};
Physical Line("outlet") = {15};
Physical Line("wall") = {10, 11, 12, 13, 14, 16};

Physical Surface("fluid") = {20};

Mesh.Smoothing = 3;
