cmake_minimum_required(VERSION 3.20)
project(NTN_OFDM_Simulator
    VERSION 1.0.0
    DESCRIPTION "C++ NTN Signal Processing and OFDM Link Simulator"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings and options
if(MSVC)
    add_compile_options(/W4 /EHsc /permissive- /utf-8)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Core library
add_library(ntn_core
    src/signal/dsp_utils.cpp
    src/signal/fft.cpp
    src/channel/noise.cpp
    src/channel/link_budget.cpp
    src/modem/modulation.cpp
    src/metrics/metrics.cpp
    src/ofdm/ofdm_modem.cpp
    src/simulation/awgn_sim.cpp
    src/synchronization/cfo.cpp
    src/channel/doppler.cpp
    src/compensation/ici_mitigation.cpp
)

target_include_directories(ntn_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Enable testing
enable_testing()

# Unit tests executable: DSP
add_executable(test_dsp tests/test_dsp.cpp)
target_link_libraries(test_dsp PRIVATE ntn_core)
add_test(NAME DSP_Unit_Tests COMMAND test_dsp)

# Unit tests executable: Modulation
add_executable(test_modulation tests/test_modulation.cpp)
target_link_libraries(test_modulation PRIVATE ntn_core)
add_test(NAME Modulation_Unit_Tests COMMAND test_modulation)

# Unit tests executable: OFDM Modem
add_executable(test_ofdm tests/test_ofdm.cpp)
target_link_libraries(test_ofdm PRIVATE ntn_core)
add_test(NAME OFDM_Unit_Tests COMMAND test_ofdm)

# Unit tests executable: AWGN Simulator
add_executable(test_awgn_sim tests/test_awgn_sim.cpp)
target_link_libraries(test_awgn_sim PRIVATE ntn_core)
add_test(NAME AWGN_Sim_Unit_Tests COMMAND test_awgn_sim)

# Unit tests executable: CFO Synchronization
add_executable(test_cfo tests/test_cfo.cpp)
target_link_libraries(test_cfo PRIVATE ntn_core)
add_test(NAME CFO_Unit_Tests COMMAND test_cfo)

# Unit tests executable: Doppler Tracking & Orbit
add_executable(test_doppler tests/test_doppler.cpp)
target_link_libraries(test_doppler PRIVATE ntn_core)
add_test(NAME Doppler_Unit_Tests COMMAND test_doppler)

# Unit tests executable: ICI Modeling & Equalization
add_executable(test_ici tests/test_ici.cpp)
target_link_libraries(test_ici PRIVATE ntn_core)
add_test(NAME ICI_Unit_Tests COMMAND test_ici)

# Example executable: Satellite Link Budget Demo
add_executable(basic_link_budget examples/basic_link_budget.cpp)
target_link_libraries(basic_link_budget PRIVATE ntn_core)

# Experiment executable: Monte Carlo AWGN Simulation Runner
add_executable(run_awgn_simulation experiments/run_awgn_simulation.cpp)
target_link_libraries(run_awgn_simulation PRIVATE ntn_core)

# Experiment executable: Carrier Frequency Offset (CFO) Experiment Runner
add_executable(run_cfo_experiment experiments/run_cfo_experiment.cpp)
target_link_libraries(run_cfo_experiment PRIVATE ntn_core)

# Experiment executable: LEO Doppler Experiment Runner
add_executable(run_doppler_experiment experiments/run_doppler_experiment.cpp)
target_link_libraries(run_doppler_experiment PRIVATE ntn_core)

# Experiment executable: Inter-Carrier Interference (ICI) Experiment Runner
add_executable(run_ici_experiment experiments/run_ici_experiment.cpp)
target_link_libraries(run_ici_experiment PRIVATE ntn_core)

# Experiment executable: Comprehensive NTN Scenarios Benchmark Matrix
add_executable(run_scenario_matrix experiments/run_scenario_matrix.cpp)
target_link_libraries(run_scenario_matrix PRIVATE ntn_core)


