// Gmsh project created on Sat Apr  8 16:13:07 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {-1, -1, 0};
Point(2) = {1, -1, 0};
Point(3) = {1, 1, 0};
Point(4) = {-1, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//+
Point(5) = {-0.5, -0.5, 0, 1.0};
Point(6) = {0.5, -0.5, 0, 1.0};
Point(7) = {-0.5, 0.5, 0, 1.0};


//+
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7};

Plane Surface(1) = {1, 2};
Plane Surface(2) = {2};


Physical Curve("bc_outer", 5) = {1, 2, 3, 4};
//Physical Curve("bc_inside", 2) = {5, 6, 7};
Physical Point("bc_outer", 2) = {5, 6, 7};

Physical Surface("mat_1", 6) = {1};
Physical Surface("mat_2", 2) = {2};

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = .6;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;


