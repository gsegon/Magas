
# Magas
add_test(NAME test_magas_2_conductors COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/2_conductors/2_conductors.json)
add_test(NAME test_magas_6_conductors_3_circuits COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/6_conductors_3_circuits/6_conductors_3_circuits.json)
add_test(NAME test_magas_EI_core COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/EI_core/EI_core.json)
add_test(NAME test_magas_BlockMagnet COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/magnet/BlockMagnet.json)
add_test(NAME test_magas_ipm_2 COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/ipm_2/ipm_2.json)
add_test(NAME test_magas_motoric_section COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/ipm_2_section/ipm_2_section.json)
add_test(NAME test_magas_motoric_section_airgap COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/ipm_2_section_airgap_band/ipm_2_section_airgap_band.json)
add_test(NAME test_magas_torque_benchmark_kelvin_1 COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/torque_benchmark_kelvin_1/torque_benchmark_kelvin_1.json)
add_test(NAME test_magas_unit_square COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/unit_square/unit_square.json)

# Magas - nonlinear
add_test(NAME test_magas_EI_core_nonlinear COMMAND magas ${CMAKE_CURRENT_SOURCE_DIR}/EI_core/EI_core_nonlinear.json)

