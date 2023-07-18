// Gmsh project created on Thu Jun 15 08:55:07 2023
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 1e-2, 0, 2*Pi};
Circle(2) = {0, 0, 0, 0.78e-2, 0, 2*Pi};
Circle(4) = {2.25e-2, 0, 0, 1.0e-2, 0, 2*Pi};
Rectangle(1) = {-0.5e-2, -0.25e-2, 0, 1e-2, 0.5e-2, 0};
Circle(9) = {2.25e-2, 0, 0, 0.001e-2, 0, 2*Pi};

Physical Curve("per_left") = {1};
Physical Curve("per_right") = {4};
Physical Surface("magnet") = {1};

Curve Loop(3) = {7, 8, 5, 6};
Curve Loop(4) = {2};
Curve Loop(5) = {7, 8, 5, 6};
Plane Surface(2) = {4, 5};
Curve Loop(6) = {1};
Curve Loop(7) = {2};
Plane Surface(3) = {6, 7};
Curve Loop(8) = {4};
Curve Loop(9) = {9};
Plane Surface(4) = {8, 9};

Transfinite Curve {1} = 20 Using Progression 1;
Periodic Curve {4} = {1};

Physical Surface("air_1") = {2};
Physical Surface("air_2") = {3};
Physical Surface("air_4") = {4};

Physical Curve("A=0") = {9};

Transfinite Curve {2} = 20 Using Progression 1;
