cmake_minimum_required(VERSION 3.16)
project(mfem_benchmarks LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if (APPLE)
  set(LIBCEED_FLAGS "-O3 -mcpu=apple-m1 -fvectorize -fslp-vectorize -ffp-contract=fast")
else()
  set(LIBCEED_FLAGS "-O3 -ffp-contract=fast")
endif()

if (ENABLE_CUDA)
  enable_language(CUDA)
  set(MFEM_USE_CUDA ON CACHE BOOL "")
  set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "")
else()
  set(MFEM_USE_CUDA OFF CACHE BOOL "")
  set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "")
endif()

if (MEFM_USE_CEED)
  set(CEED_DIR ${PROJECT_SOURCE_DIR}/libCEED CACHE STRING "")
endif()

add_subdirectory(mfem)
#add_dependencies(mfem libCEED)

if (ENABLE_ASAN)
  add_link_options(-fsanitize=address)
  add_compile_options(-fsanitize=address)
endif()

file(GLOB_RECURSE cpp_tests ${PROJECT_SOURCE_DIR}/tests/*.cpp)

# note: this is an unnatural way to handle CUDA compilation, 
# but the way mfem is written sort of forces us to do it
if (ENABLE_CUDA)
  set_source_files_properties(${cpp_tests} PROPERTIES LANGUAGE CUDA)
  message("${cpp_tests}")
endif()

foreach(filename ${cpp_tests})		

  get_filename_component(testname ${filename} NAME_WE)

  add_executable(${testname} ${filename})
  target_link_libraries(${testname} PUBLIC mfem)

  add_test(${testname} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${testname})		

endforeach(filename ${cpp_tests})