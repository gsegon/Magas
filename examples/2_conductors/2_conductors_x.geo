// Gmsh project created on Fri Apr 21 16:21:20 2023
SetFactory("OpenCASCADE");


//+
Circle(5) = {-0.25, 0, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {0.25, 0.0, 0, 0.1, 0, 2*Pi};
//+
//+

//+
Point(7) = {-1, -1, 0, 1.0};
//+
Point(8) = {1, -1, 0, 1.0};
//+
Point(9) = {1, 1, 0, 1.0};
//+
Point(10) = {-1, 1, 0, 1.0};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 7};
//+
Line Loop(1) = {5};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {6};
//+
Plane Surface(2) = {2};
//
Line Loop(3) = {7, 8, 9, 10};
Plane Surface(3) = {3, 1, 2};

Physical Line("dc_boundary", 100) = {7, 8, 9, 10};

Physical Surface("conductor_1", 1) = {1};
Physical Surface("conductor_2", 2) = {2};
Physical Surface("air", 3) = {3};


//+
//Physical Line("bc_conductor_1", 1) = {5};
