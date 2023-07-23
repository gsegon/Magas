// Gmsh project created on Thu May  4 12:33:04 2023
SetFactory("OpenCASCADE");
ShapeFromFile("motoric_section_airgap.brep")

//Line(241) = {187, 5};
//Line(242) = {186, 1};

//+
Point(340) = {46.2, 0, 0, 1.0};
Point(341) = {0, 46.2, 0, 1.0};
//+
Point(342) = {46.4, 0, 0, 1.0};
Point(343) = {0, 46.4, 0, 1.0};

Point(344) = {0, 0, 0, 1.0};
//+
Circle(419) = {340, 344, 341};
//+
Circle(420) = {342, 344, 343};


//+
Line(421) = {186, 340};
//+
Line(422) = {340, 342};
//+
Line(423) = {342, 1};
//+
Line(424) = {187, 341};
//+
Line(425) = {341, 343};
//+
Line(426) = {343, 5};

//+
Curve Loop(23) = {186, 424, -419, -421};
//+
Plane Surface(19) = {23};
//+
Curve Loop(24) = {425, -420, -422, 419};
//+
Plane Surface(20) = {24};

Coherence;//+
Curve Loop(25) = {520, -519, 518, -531, -532, 512, -511, 510, -535, -536, 504, -503, 502, -539, -540, 496, -495, 494, -543, -544, 488, -487, 486, -547, -548, 480, -479, 478, -533, -534, 472, -471, 470, -529, -530, 464, -463, 462, -541, -542, 456, -455, 454, -537, -538, 448, -447, 446, -549, -550, 440, -439, 438, -545, -546, 432, 431, -426, -420, 423, 427, 526, -527, -528};
//+
Plane Surface(21) = {25};
//+
Curve Loop(27) = {590, -589, -599, 593, 591};
//+
Plane Surface(22) = {27};
//+
Curve Loop(29) = {598, -584, 585};
//+
Plane Surface(23) = {29};
//+
Curve Loop(30) = {574, -600, -583};
//+
Plane Surface(24) = {30};
//+
Curve Loop(31) = {601, 579, -578, -577, -576};
//+
Plane Surface(25) = {31};
//+
Curve Loop(32) = {560, -559, -595, 563, 561};
//+
Plane Surface(26) = {32};
//+
Curve Loop(34) = {594, -554, 555};
//+
Plane Surface(27) = {34};
//+
Curve Loop(35) = {565, -596, -564};
//+
Plane Surface(28) = {35};
//+
Curve Loop(36) = {597, 570, -569, -568, -567};
//+
Plane Surface(29) = {36};
//+
Physical Curve("A0", 505) = {553, 429};

Periodic Curve {551} = {552};
Periodic Curve {421} = {424};
Periodic Curve {422} = {425};
Periodic Curve {423} = {426};
Periodic Curve {428} = {430};

Physical Curve("bc_periodic_a", 518) = {551, 421, 422, 423, 428};
Physical Curve("bc_periodic_b", 519) = {552, 424, 425, 426, 430};

Physical Surface("Magnet_N_1", 512) = {17};
Physical Surface("Magnet_N_2", 513) = {18};
Physical Surface("Magnet_S_1", 514) = {15};
Physical Surface("Magnet_S_2", 515) = {16};

Physical Surface("Coil_A+", 506) = {11, 2};
Physical Surface("Coil_A-", 507) = {5, 12};
Physical Surface("Coil_B-", 508) = {13, 7};
Physical Surface("Coil_C+", 509) = {3, 9};
Physical Surface("Coil_B+", 510) = {8, 10};
Physical Surface("Coil_C-", 511) = {6, 4};

Physical Surface("AirPockets", 516) = {22, 23, 24, 25, 26, 27, 28, 29};
Physical Surface("Air", 2000) = {19};
Physical Surface("Air2", 2001) = {21};
Physical Surface("AirGap", 517) = {20};//+
//+
Physical Surface("RotorCore", 601) = {14};
//+
Physical Surface("StatorCore", 602) = {1};
//+
Transfinite Curve {420} = 90 Using Progression 1;
//+
Transfinite Curve {186} = 30 Using Progression 1;
//+
Transfinite Curve {419} = 90 Using Progression 1;
//+