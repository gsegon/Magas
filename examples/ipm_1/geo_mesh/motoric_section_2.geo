// Gmsh project created on Thu May  4 12:33:04 2023
SetFactory("OpenCASCADE");
ShapeFromFile("ipm_export_section.brep")

//+
Line(164) = {149, 95};
//+
Line(165) = {152, 65};
//+
//+
Coherence;
//+
Curve Loop(1) = {281, 266, 267};
//+
Plane Surface(20) = {1};
//+
Curve Loop(2) = {263, 264, 280};
//+
Plane Surface(21) = {2};
//+
Curve Loop(3) = {272, -282, 271};
//+
Plane Surface(22) = {3};
//+
Curve Loop(4) = {283, -275, -274};
//+
Plane Surface(23) = {4};
//+
Curve Loop(6) = {164, -257, 171, 172, 173, 174, 175, -255, 183, 184, 185, 186, 187, -253, 195, 196, 197, 198, 199, -251, 207, 208, 209, 210, 211, -249, 219, 220, 221, 222, 223, -247, 231, 232, 233, 234, 235, -244, -165, 279};
//+
Plane Surface(24) = {6};
//+
Physical Curve("a=0", 505) = {241};
//+
Physical Surface("air_1", 200) = {19};
//+
Physical Surface("air_pockets", 201) = {20, 22, 21, 23};
//+
Physical Surface("rotor_core", 3) = {16};
//+
Physical Surface("stator_core", 4) = {1};
//+
Physical Surface("phase_a+", 10) = {8, 7, 13, 12};
//+
Physical Surface("phase_b-", 11) = {6, 5, 11, 10};
//+
Physical Surface("phase_c+", 12) = {4, 2, 3};
//+
Physical Surface("phase_c-", 13) = {9, 15, 14};
//+
Physical Surface("airgap") = {24};

Physical Curve("per_1", 507) = {166, 259, 256, 164, 276, 285};
//+
Physical Curve("per_2", 508) = {240, 245, 242, 165, 278, 284};

Periodic Line {166} = {240};
Periodic Line {259} = {245};
Periodic Line {256} = {242};
Periodic Line {164} = {165};
Periodic Line {276} = {278};
Periodic Line {285} = {284};

Physical Surface("magnet_1") = {17};
//+
Physical Surface("magnet_2") = {18};

Transfinite Curve {279} = 90 Using Progression 1;
Mesh 2;//+

