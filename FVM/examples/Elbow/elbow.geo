// Structured elbow geometry from the icoFoam tutorial
//
//               -Ly- 
//               .---. -
//               |(3)| L
//               |   | x
//    --Lx--    .´---| -   ^ y 
//  -  _____..·´ (2).´     |
//  L | (1) :    .  |      |
//  y |_____:.··´|__|      ·-----> 
//  -            -lx-            x
//

Ly = 12; // main pipe diameter
Lx = 24; // main pipe length

lx = 4; // secondary pipe diameter

// Dimensions
ndim_x = 25; // x-direction (sections 1 and 3)
ndim_x2 = 18; // x-direction in elbow section (section 2)
ndim_y = 25; // y-direction (all sections)
ndim_pipe = 15; // secondary pipe direction

// progression and bump increase ratio
prog = 1.1; // (sections 1 and 3)
prog2 = 1.2; // elbow section (2)
bump = 0.25; // all sections

// Calculations
Lx2 = Lx + Ly;
Ly2 = 2*Ly;
angle = Asin((lx/2)*Sqrt(2)/(Lx));

// Section (1)
Point(1) = {0, 0, 0, 1};
Point(2) = {Lx, 0, 0, 1};
Point(3) = {Lx, Ly, 0, 1};
Point(4) = {0, Ly, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Section (3)
Point(5) = {Lx2, Ly2, 0, 1};
Point(6) = {Lx2, Ly2+Lx, 0, 1};
Point(7) = {Lx2+Ly, Ly2+Lx, 0, 1};
Point(8) = {Lx2+Ly, Ly2, 0, 1};

Line(5) = {6, 5};
Line(6) = {7, 6};
Line(7) = {8, 7};
Line(8) = {5, 8};

// Section (2) -> split into four surfaces 
Point(9) = {Lx, 2*Ly, 0, 1};
Point(10) = {Lx+Ly*Cos(Pi/4-angle), Lx-Ly*Sin(Pi/4-angle), 0, 1};
Point(11) = {Lx+Ly*Cos(Pi/4+angle), Lx-Ly*Sin(Pi/4+angle), 0, 1};
Point(12) = {Lx+Lx*Cos(Pi/4-angle), Lx-Lx*Sin(Pi/4-angle), 0, 1};
Point(13) = {Lx+Lx*Cos(Pi/4+angle), Lx-Lx*Sin(Pi/4+angle), 0, 1};

Circle(9) = {3, 9, 11};
Circle(10) = {10, 9, 11};
Circle(11) = {10, 9, 5};
Circle(12) = {13, 9, 2};
Circle(13) = {8, 9, 12};

Point(14) = {Lx+Lx*Cos(Pi/4+angle),0 , 0, 1};
Point(15) = {Lx+Lx*Cos(Pi/4-angle),0 , 0, 1};

Line(14) = {13, 12};
Line(15) = {14, 13};
Line(16) = {15, 14};
Line(17) = {12, 15};

Line(18) = {12, 10};
Line(19) = {11, 13};

// Create surfaces -- 2D
Line Loop(20) = {1, 2, 3, 4};
Plane Surface(21) = {20};

Line Loop(22) = {5, 8, 7, 6};
Plane Surface(23) = {22};

Line Loop(24) = {15, 14, 17, 16};
Plane Surface(25) = {24};

Line Loop(26) = {18, 10, 19, 14};
Plane Surface(27) = {26};

Line Loop(28) = {2, 9, 19, 12};
Plane Surface(29) = {28};

Line Loop(30) = {18, 13, 8, 11};
Plane Surface(31) = {30};

Transfinite Surface {21, 23, 25, 27, 29, 31};
Recombine Surface {21, 23, 25, 27, 29, 31};

Transfinite Line {1, 5} = ndim_x + 1 Using Progression prog;
Transfinite Line {3, 7} = ndim_x + 1 Using Progression 1./prog;
Transfinite Line {11, 12} = ndim_x2 + 1 Using Progression prog2;
Transfinite Line {9, 13} = ndim_x2 + 1 Using Progression 1./prog2;
Transfinite Line {10, 14, 16} = ndim_pipe + 1 Using Bump bump;
Transfinite Line {15, 17} = ndim_x2 + 1 Using Bump bump;
Transfinite Line {2, 4, 6, 8, 18, 19} = ndim_y + 1 Using Bump bump;

// Define boundary surfaces
Physical Line("wall") = {1, 3, 9, 10, 11, 5, 7, 13, 17, 15, 12, 1};
Physical Line("outlet") = {6};
Physical Line("main-inlet") = {4};
Physical Line("secondary-inlet") = {16};

Physical Surface("fluid") = {21, 23, 25, 27, 29, 31};
