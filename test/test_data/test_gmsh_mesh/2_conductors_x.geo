// Gmsh project created on Fri Apr 21 16:21:20 2023
SetFactory("OpenCASCADE");


//+
Circle(1) = {-0.25, 0, 0, 0.1, 0, 2*Pi};
//+
Circle(2) = {0.25, 0.0, 0, 0.1, 0, 2*Pi};
//+
//+

//+
Point(3) = {-1, -1, 0, 1.0};
//+
Point(4) = {1, -1, 0, 1.0};
//+
Point(5) = {1, 1, 0, 1.0};
//+
Point(6) = {-1, 1, 0, 1.0};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 3};
//+
Line Loop(1) = {1};
//+
Plane Surface(201) = {1};
//+
Line Loop(2) = {2};
//+
Plane Surface(202) = {2};
//


//Physical Line("dc_boundary", 100) = {3, 4, 5, 6};

//Physical Surface("conductor_1", 1) = {201};
//Physical Surface("conductor_2", 2) = {202};
//Physical Surface("air", 3) = {203};


//+
//Physical Curve("dc_conductor_1", 101) = {1};
//Physical Point("dc_conductor_1", 101) = {1};

//Physical Curve("dc_conductor_2", 102) = {2};
//Physical Point("dc_conductor_2", 102) = {2};
//+
Curve Loop(3) = {6, 3, 4, 5};
//+
Curve Loop(4) = {6, 3, 4, 5};
//+
Plane Surface(203) = {4};
//+
BooleanDifference{ Surface{203}; Delete; }{ Surface{201}; Surface{202}; }

