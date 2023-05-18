// Gmsh project created on Sun Apr 23 11:40:08 2023
SetFactory("OpenCASCADE");

w_bound = 2*1e-3;
h_bound = 2*1e-3;

w_magnet = 0.25*1e-3;
h_magnet = 0.5*1e-3;

//+
Point(1) = {-w_bound, -h_bound, 0, 1.0};
//+
Point(2) = {w_bound, -h_bound, 0, 1.0};
//+
Point(3) = {w_bound, h_bound, 0, 1.0};
//+
Point(4) = {-w_bound, h_bound, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Point(5) = {-w_magnet, -h_magnet, 0, 1.0};
//+
Point(6) = {w_magnet, -h_magnet, 0, 1.0};
//+
Point(7) = {w_magnet, h_magnet, 0, 1.0};
//+
Point(8) = {-w_magnet, h_magnet, 0, 1.0};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {8, 5, 6, 7};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {8, 5, 6, 7};
//+
Plane Surface(2) = {3};

Physical Surface("Air", 1) = {1};
Physical Surface("Magnet", 2) = {2};
Physical Curve("A0", 505) = {1, 2, 3, 4};
