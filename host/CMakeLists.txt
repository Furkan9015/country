cmake_minimum_required(VERSION 3.10)   # safer than 3.5 now
project(upvc)

if (NOT DEFINED UPMEM_HOME)
        if ( "$ENV{UPMEM_HOME}" STREQUAL "")
                set(UPMEM_HOME "/usr")
        else ()
                set(UPMEM_HOME $ENV{UPMEM_HOME})
        endif ()
endif ()

include(${UPMEM_HOME}/share/upmem/cmake/include/host/DpuHost.cmake)

# Find OpenMP for parallel optimizations
find_package(OpenMP)

set(DPU_PROJECT_RELATIVE_PATH ../dpu)
set(DPU_BINARY_NAME dpu_task)
set(NR_TASKLETS 16)
set(CMAKE_C_FLAGS "--std=gnu99 -O3 -march=native -flto -Wall -Wextra  -g3 -DNR_TASKLETS=${NR_TASKLETS} -DDPU_BINARY=\\\"${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}/${DPU_BINARY_NAME}\\\"")

# Enable LTO for linker
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")

# Add OpenMP flags if available
if(OPENMP_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
else()
    message(WARNING "OpenMP not found - parallel optimizations will be disabled")
endif()
link_directories("${DPU_HOST_LINK_DIRECTORIES}")

# WFA2 library paths
set(WFA2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external/WFA2-lib)
set(WFA2_INCLUDE ${WFA2_DIR})
set(WFA2_LIB ${WFA2_DIR}/lib/libwfa.a)

file(GLOB_RECURSE SOURCES src/*.c)

# Traditional UPVC build (baseline)
add_executable(upvc ${SOURCES})
target_include_directories(upvc PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/ ${WFA2_INCLUDE})
target_link_libraries(upvc ${DPU_HOST_LIBRARIES} pthread ${WFA2_LIB} m)

# OCOCO-enabled UPVC build (Stage 1)
add_executable(upvc_ococo ${SOURCES})
target_include_directories(upvc_ococo PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/ ${WFA2_INCLUDE})
target_compile_definitions(upvc_ococo PRIVATE USE_OCOCO_COUNTERS)
target_link_libraries(upvc_ococo ${DPU_HOST_LIBRARIES} pthread ${WFA2_LIB} m)

# WFA2-accelerated UPVC build (Option A)
add_executable(upvc_wfa2 ${SOURCES})
target_include_directories(upvc_wfa2 PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/ ${WFA2_INCLUDE})
target_compile_definitions(upvc_wfa2 PRIVATE USE_WFA2)
target_link_libraries(upvc_wfa2 ${DPU_HOST_LIBRARIES} pthread ${WFA2_LIB} m)

# Stage 2: OCOCO Lite (counters-only, no variant hash tables)
add_executable(upvc_stage2 ${SOURCES})
target_include_directories(upvc_stage2 PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/ ${WFA2_INCLUDE})
target_compile_definitions(upvc_stage2 PRIVATE USE_OCOCO_COUNTERS OCOCO_LITE_MODE)
target_link_libraries(upvc_stage2 ${DPU_HOST_LIBRARIES} pthread ${WFA2_LIB} m)

# Stage 2.5: OCOCO Hybrid (hash tables + counter evidence at variant positions)
add_executable(upvc_hybrid ${SOURCES})
target_include_directories(upvc_hybrid PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/ ${WFA2_INCLUDE})
target_compile_definitions(upvc_hybrid PRIVATE USE_OCOCO_COUNTERS OCOCO_HYBRID_MODE)
target_link_libraries(upvc_hybrid ${DPU_HOST_LIBRARIES} pthread ${WFA2_LIB} m)

# Stage 2.6: OCOCO Bayesian Calling (hybrid mode + Bayesian genotype calling)
add_executable(upvc_bayesian ${SOURCES})
target_include_directories(upvc_bayesian PUBLIC "${DPU_HOST_INCLUDE_DIRECTORIES}" inc/ ../common/inc/ ${WFA2_INCLUDE})
target_compile_definitions(upvc_bayesian PRIVATE USE_OCOCO_COUNTERS OCOCO_HYBRID_MODE OCOCO_BAYESIAN_CALLING)
target_link_libraries(upvc_bayesian ${DPU_HOST_LIBRARIES} pthread ${WFA2_LIB} m)

set(NB_DPU_MARK)
if (NB_DPU)
        set(NB_DPU_MARK "-n")
endif()

## Traditional UPVC test targets
add_custom_target(mrams yes | ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g index -n 128
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(run yes | ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map ${NB_DPU_MARK} ${NB_DPU} -f
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(run_simu yes | ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map -f -s
        DEPENDS upvc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(check ${CMAKE_CURRENT_SOURCE_DIR}/../tests/compareVCF.py ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_upvc_ref.vcf ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_upvc.vcf -c
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

## OCOCO test targets (Stage 1)
add_custom_target(run_ococo yes | ${CMAKE_CURRENT_BINARY_DIR}/upvc_ococo -i chr22 -g map ${NB_DPU_MARK} ${NB_DPU} -f
        DEPENDS upvc_ococo
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(run_ococo_simu yes | ${CMAKE_CURRENT_BINARY_DIR}/upvc_ococo -i chr22 -g map -f -s
        DEPENDS upvc_ococo
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)
add_custom_target(check_ococo ${CMAKE_CURRENT_SOURCE_DIR}/../tests/compareVCF.py ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_upvc_ref.vcf ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration/chr22_upvc_ococo.vcf -c
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

## Benchmark: Compare OCOCO vs Baseline
add_custom_target(benchmark_stage1
        COMMAND echo "=== Running Baseline UPVC ==="
        COMMAND /usr/bin/time -v ${CMAKE_CURRENT_BINARY_DIR}/upvc -i chr22 -g map -f -s > baseline.log 2>&1
        COMMAND echo "=== Running OCOCO UPVC (Stage 1) ==="
        COMMAND /usr/bin/time -v ${CMAKE_CURRENT_BINARY_DIR}/upvc_ococo -i chr22 -g map -f -s > ococo_stage1.log 2>&1
        COMMAND echo "=== Comparing Results ==="
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../tests/compareVCF.py chr22_upvc_ref.vcf chr22_upvc.vcf -c
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../tests/compareVCF.py chr22_upvc_ref.vcf chr22_upvc_ococo.vcf -c
        COMMAND echo "=== Performance Summary ==="
        COMMAND grep "Maximum resident" baseline.log ococo_stage1.log || true
        COMMAND grep "Elapsed" baseline.log ococo_stage1.log || true
        DEPENDS upvc upvc_ococo
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../tests/chr22_integration)

include(ExternalProject)
ExternalProject_Add(
        ${DPU_BINARY_NAME}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${DPU_PROJECT_RELATIVE_PATH}
        CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE=${UPMEM_HOME}/share/upmem/cmake/dpu.cmake -DUPMEM_HOME=${UPMEM_HOME} -DNR_TASKLETS=${NR_TASKLETS}
        BUILD_ALWAYS TRUE
        INSTALL_COMMAND ""
)

add_dependencies(upvc dpu_task)
add_dependencies(upvc_ococo dpu_task)
add_dependencies(upvc_stage2 dpu_task)
add_dependencies(upvc_hybrid dpu_task)
add_dependencies(upvc_bayesian dpu_task)
