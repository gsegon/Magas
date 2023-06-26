// Gmsh project created on Sat Apr  8 16:13:07 2023
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-1.0,-1.0, 0, 2, 2, 0};
//+
Physical Curve("dc_0", 5) = {4, 3, 2, 1};
//+
Physical Surface("domain", 6) = {1};
//+

