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
Rectangle(1) = {-0.5e-2, -0.25e-2, 0, 1e-2, 0.5e-2, 0};
//+

//+
//+
Physical Curve("per_left") = {1};
//+
Physical Curve("per_right") = {4};
//+
Physical Surface("magnet") = {1};
//+
Curve Loop(2) = {3};
//+
Curve Loop(3) = {7, 8, 5, 6};
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

//+
Physical Surface("air_1") = {2};

Physical Surface("air_2") = {3};

Physical Surface("air_3") = {4};


Periodic Curve {2} = {3};


//+
Circle(9) = {2.25e-2, 0, 0, 0.001e-2, 0, 2*Pi};
//+
Curve Loop(9) = {4};
//+
Curve Loop(10) = {9};
//+
Plane Surface(5) = {9, 10};
//+
Physical Curve("A=0") = {9};
//+
Physical Surface("air_4") = {5};
