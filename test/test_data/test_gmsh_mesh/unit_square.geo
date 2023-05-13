// Gmsh project created on Sat Apr  8 16:13:07 2023
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-1.0,-1.0, 0, 2, 2, 0};

//+
Point(5) = {-0.5, 0.5, 0, 1.0};
//+
Point(6) = {0.5, -0.5, 0, 1.0};
//+
Line(5) = {5, 6};

Line{5} In Surface{1};

Physical Curve("dc_0", 5) = {4, 3, 2, 1};
Physical Surface("domain", 6) = {1};
Physical Curve("bc_line", 2) = {5};

