// Gmsh project created on Thu Jun 15 08:55:07 2023
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 1e-2, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 0.78e-2, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 0.72e-2, 0, 2*Pi};
//+
Circle(4) = {2.25e-2, 0, 0, 1.0e-2, 0, 2*Pi};
//+
Circle(5) = {0, 0, 0, 0.25e-2, 0, 2*Pi};
//+
Circle(9) = {2.25e-2, 0, 0, 0.001e-2, 0, 2*Pi};
//+
Curve Loop(1) = {5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3};
//+
Curve Loop(3) = {5};
//+
Plane Surface(2) = {2, 3};
//+
Curve Loop(4) = {2};
//+
Curve Loop(5) = {3};
//+
Plane Surface(3) = {4, 5};
//+
Curve Loop(6) = {1};
//+
Curve Loop(7) = {2};
//+
Plane Surface(4) = {6, 7};
//+
Curve Loop(8) = {4};
//+
Curve Loop(9) = {9};
//+
Plane Surface(5) = {8, 9};

Physical Curve("per_left", 1) = {1};
Physical Curve("per_right", 2) = {4};
Physical Surface("air_1", 4) = {2};
Physical Surface("air_2", 5) = {3};
Physical Surface("air_3", 6) = {4};
Physical Curve("A=0", 7) = {9};
Physical Surface("air_4", 8) = {5};
Physical Surface("Conductor", 3) = {1};
//+
Transfinite Curve {5} = 100 Using Progression 1;
