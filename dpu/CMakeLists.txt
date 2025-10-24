cmake_minimum_required(VERSION 3.10)   # safer than 3.5 now

project(dpu C ASM)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)

# Option to enable Stage 2 DPU counter updates
# Temporarily hardcoded to ON for testing
set(CMAKE_C_FLAGS "-O2 -g -fstack-size-section -DNR_TASKLETS=${NR_TASKLETS} -DSTACK_SIZE_DEFAULT=192 -DUSE_DPU_COUNTERS")
message(STATUS "Stage 2: DPU counter updates ENABLED")

INCLUDE_DIRECTORIES(inc)
INCLUDE_DIRECTORIES(../common/inc/)

set(SOURCES_OPT2 src/result_pool.c src/request_pool.c src/dout.c src/nodp.c src/odpd.c src/odpd_init_opt.c src/odpd_opt.S src/task.c src/nodp.S src/edit_distance.c src/mram_counters.c)

add_executable(dpu_task ${SOURCES_OPT2})
