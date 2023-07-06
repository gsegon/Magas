
# Magas
add_test(NAME test_magas_2_conductors COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/2_conductors/2_conductors.json)
add_test(NAME test_magas_6_conductors_3_circuits COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/6_conductors_3_circuits/6_conductors_3_circuits.json)
add_test(NAME test_magas_EI_core COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/EI_core/EI_core.json)
add_test(NAME test_magas_BlockMagnet COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/magnet/BlockMagnet.json)
add_test(NAME test_magas_motoric COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/motoric/motoric.json)
add_test(NAME test_magas_motoric_section COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/motoric_section/motoric_section.json)
add_test(NAME test_magas_motoric_section_airgap COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/motoric_section_airgap/motoric_section_airgap.json)
add_test(NAME test_magas_torque_benchmark_kelvin_1 COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.json)
add_test(NAME test_magas_unit_square COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/unit_square/unit_square.json)

# amnl
add_test(NAME test_amnl_unit_square COMMAND amnl ${CMAKE_CURRENT_SOURCE_DIR}/unit_square/unit_square.json)
add_test(NAME test_amnl_EI_core COMMAND amnl ${CMAKE_CURRENT_SOURCE_DIR}/EI_core/EI_core.json)
add_test(NAME test_amnl_EI_core_nonlinear COMMAND amnl ${CMAKE_CURRENT_SOURCE_DIR}/EI_core/EI_core_nonlinear.json)

