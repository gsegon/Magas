// Gmsh project created on Fri Jun  9 12:58:15 2023
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-2, -2, 0, 4, 4, 0};
//+
Circle(5) = {-0.5, 1, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {0.5, 1.2, 0, 0.1, 0, 2*Pi};
//+
Circle(7) = {0.8, 0.6, 0, 0.1, 0, 2*Pi};
//+
Circle(8) = {0.5, -0.5, 0, 0.1, 0, 2*Pi};
//+
Circle(9) = {-0.3, -0.6, 0, 0.1, 0, 2*Pi};
//+
Circle(10) = {-0.8, -0.2, 0, 0.1, 0, 2*Pi};
//+
Curve Loop(2) = {5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {9};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {10};
//+
Plane Surface(7) = {7};
//+

//+
//+
Coherence;
//+
Physical Curve("A=0") = {2, 1, 3, 4};
//+
Physical Surface("air") = {8};
//+
Physical Surface("J2in") = {2};
//+
Physical Surface("J1out") = {3};
//+
Physical Surface("J3out") = {4};
//+
Physical Surface("J3in") = {5};
//+
Physical Surface("J2out") = {6};
//+
Physical Surface("J1in") = {7};
