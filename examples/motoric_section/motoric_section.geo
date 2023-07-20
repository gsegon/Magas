// Gmsh project created on Thu May  4 12:33:04 2023
SetFactory("OpenCASCADE");
ShapeFromFile("motoric_section.brep")

Line(241) = {187, 5};
Line(242) = {186, 1};

Coherence;

//+
Curve Loop(1) = {408, 407, -406, -416, 410};
//+
Plane Surface(19) = {1};
//+
Curve Loop(3) = {415, -401, 402};
//+
Plane Surface(20) = {3};
//+
Curve Loop(4) = {391, -417, -400};
//+
Plane Surface(21) = {4};
//+
Curve Loop(5) = {418, 396, -395, -394, -393};
//+
Plane Surface(22) = {5};
//+
Curve Loop(6) = {377, -376, -412, 380, 378};
//+
Plane Surface(23) = {6};
//+
Curve Loop(8) = {411, -371, 372};
//+
Plane Surface(24) = {8};
//+
Curve Loop(9) = {382, -413, -381};
//+
Plane Surface(25) = {9};
//+
Curve Loop(10) = {414, 387, -386, -385, -384};
//+
Plane Surface(26) = {10};
//+
Curve Loop(11) = {241, -247, -248, 362, 361, -254, 255, -256, 366, 365, -262, 263, -264, 354, 353, -270, 271, -272, 358, 357, -278, 279, -280, 346, 345, -286, 287, -288, 350, 349, -294, 295, -296, 364, 363, -302, 303, -304, 360, 359, -310, 311, -312, 356, 355, -318, 319, -320, 352, 351, -326, 327, -328, 348, 347, -334, 335, -336, 344, 343, -342, -243, -242, 368};
//+
Plane Surface(27) = {11};

//Transfinite Curve {242} = 1 Using Progression 1;
//Transfinite Curve {241} = 1 Using Progression 1;
Periodic Line {246} = {244};
Periodic Line {241} = {242};
Periodic Line {369} = {367};

Transfinite Curve {368} = 60 Using Progression 1;

Mesh 2;
RecombineMesh;
Mesh.SubdivisionAlgorithm = 1;
//RefineMesh;

//

//+
//
//Physical Curve("A0", 505) = {1};
//Physical Surface("Core_Rotor", 1) = {1};
//Physical Surface("Core_Stator", 2) = {2};
//+
Physical Curve("A0", 505) = {370, 245};
//+
Physical Surface("Core_Rotor", 1) = {14};
//+
Physical Surface("Core_Stator", 2) = {1};
//+
Physical Surface("Coil_A+", 506) = {2, 11};
//+
Physical Surface("Coil_A-", 507) = {5, 12};
//+
Physical Surface("Coil_B-", 508) = {13, 7};
//+
Physical Surface("Coil_C+", 509) = {9, 3};
//+
Physical Surface("Coil_B+", 510) = {10, 8};

Physical Surface("Coil_C-", 511) = {6, 4};
//+
Physical Surface("Magnet_N_1") = {17};
//+
Physical Surface("Magnet_N_2") = {18};
//+
Physical Surface("Magnet_S_1") = {15};
//+
Physical Surface("Magnet_S_2") = {16};
//+
Physical Surface("AirPockets") = {19, 20, 21, 22, 23, 24, 25, 26};
//+
Physical Surface("Airgap") = {27};
////+
//Physical Curve("bc_periodic_rotor_a") = {369};
//Physical Curve("bc_periodic_rotor_b") = {367};
////+
//Physical Curve("bc_periodic_stator_a") = {246};
//Physical Curve("bc_periodic_stator_b") = {244};
////+
//Physical Curve("bc_periodic_airgap_a") = {241};
//Physical Curve("bc_periodic_airgap_b") = {242};

//+
Physical Curve("bc_periodic_a") = {369, 241, 246};
Physical Curve("bc_periodic_b") = {367, 242, 244};
//+

