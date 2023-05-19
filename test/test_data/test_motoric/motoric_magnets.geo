Merge "motoric.brep";
SetFactory("OpenCASCADE");
Circle(932) = {0, 0, 0, 170/2, 0, 2*Pi};
//+
Curve Loop(85) = {932};
//+
Plane Surface(67) = {85};
//+
Coherence;
//+
Physical Curve("A0", 505) = {1};
Physical Surface("Core_Rotor", 1) = {1};
Physical Surface("Core_Stator", 2) = {2};

Physical Surface("Coil_1_a+", 101) = {3, 5};
Physical Surface("Coil_1_a-", 102) = {4, 6};
Physical Surface("Coil_1_b-", 103) = {7, 9};
Physical Surface("Coil_1_b+", 104) = {8, 10};
Physical Surface("Coil_1_c+", 105) = {11, 13};
Physical Surface("Coil_1_c-", 106) = {12, 14};

Physical Surface("Coil_2_a+", 107) = {3+12, 5+12};
Physical Surface("Coil_2_a-", 108) = {4+12, 6+12};
Physical Surface("Coil_2_b-", 109) = {7+12, 9+12};
Physical Surface("Coil_2_b+", 110) = {8+12, 10+12};
Physical Surface("Coil_2_c+", 111) = {11+12, 13+12};
Physical Surface("Coil_2_c-", 112) = {12+12, 14+12};

Physical Surface("Coil_3_a+", 113) = {3+2*12, 5+2*12};
Physical Surface("Coil_3_a-", 114) = {4+2*12, 6+2*12};
Physical Surface("Coil_3_b-", 115) = {7+2*12, 9+2*12};
Physical Surface("Coil_3_b+", 116) = {8+2*12, 10+2*12};
Physical Surface("Coil_3_c+", 117) = {11+2*12, 13+2*12};
Physical Surface("Coil_3_c-", 118) = {12+2*12, 14+2*12};

Physical Surface("Coil_4_a+", 119) = {3+3*12, 5+3*12};
Physical Surface("Coil_4_a-", 120) = {4+3*12, 6+3*12};
Physical Surface("Coil_4_b-", 121) = {7+3*12, 9+3*12};
Physical Surface("Coil_4_b+", 122) = {8+3*12, 10+3*12};
Physical Surface("Coil_4_c+", 123) = {11+3*12, 13+3*12};
Physical Surface("Coil_4_c-", 124) = {12+3*12, 14+3*12};

Physical Surface("Magnet1_n", 651) = {51};
Physical Surface("Magnet2_n", 652) = {52};
Physical Surface("Magnet3_s", 653) = {53};
Physical Surface("Magnet4_s", 654) = {54};
Physical Surface("Magnet5_n", 655) = {55};
Physical Surface("Magnet6_n", 656) = {56};
Physical Surface("Magnet7_s", 657) = {57};
Physical Surface("Magnet8_s", 658) = {58};
Physical Surface("Magnet9_n", 659) = {59};
Physical Surface("Magnet10_n", 660) = {60};
Physical Surface("Magnet11_s", 661) = {61};
Physical Surface("Magnet12_s", 662) = {62};
Physical Surface("Magnet13_n", 663) = {63};
Physical Surface("Magnet14_n", 664) = {64};
Physical Surface("Magnet15_s", 665) = {65};
Physical Surface("Magnet16_s", 666) = {66};

Physical Surface("air", 127) = {69,
70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
100};

Physical Surface("air_gap", 128) = {101};
Physical Surface("Shaft", 129) = {68};
Physical Surface("AirOutside", 130) = {67};






