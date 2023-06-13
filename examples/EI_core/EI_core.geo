// Gmsh project created on Sat Apr 15 18:45:06 2023
//+

SetFactory("OpenCASCADE");

inch_to_mm = 25.4;
mm_to_m = 1e-3;

width_coil= 0.25*inch_to_mm*mm_to_m;
height_coil = 0.5*inch_to_mm*mm_to_m;
width_i_section = 1.5*inch_to_mm*mm_to_m;
height_i_section = 0.25*inch_to_mm*mm_to_m;
height_e_section = 0.75*inch_to_mm*mm_to_m;
height_airgap = 0.025*inch_to_mm*mm_to_m;

Rectangle(1) = {-width_i_section, -width_i_section, 0, 2*width_i_section, 2*width_i_section, 0};


Rectangle(2) = {-width_i_section/2.0, -height_i_section-height_airgap/2, 0, width_i_section, height_i_section, 0};
Rectangle(3) = {-width_i_section/2.0, -height_airgap/2, 0, width_i_section, height_airgap, 0};
Rectangle(4) = {-width_i_section/2.0, height_airgap/2, 0, width_i_section, height_e_section, 0};
Rectangle(5) = {-2*width_coil, height_airgap/2, 0, width_coil, height_airgap/2+height_coil, 0};
Rectangle(6) = {1*width_coil, height_airgap/2, 0, width_coil, height_airgap/2+height_coil, 0};

Coherence;


Physical Curve("a=0", 44) = {23, 25, 24, 22};
Physical Surface("e_core", 200) = {8};
Physical Surface("i_core", 201) = {2};
Physical Surface("coil_1", 202) = {5};
Physical Surface("coil_2", 203) = {6};
Physical Surface("air_surrounding", 204) = {7};
Physical Surface("air_gap", 205) = {3};
Physical Point("a=0", 44) = {25, 26, 24, 23};

Transfinite Curve {42, 17, 38, 21, 34} = 30 Using Progression 1;
Transfinite Curve {43} = 120 Using Progression 1;

