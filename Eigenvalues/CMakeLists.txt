cmake_minimum_required(VERSION 3.10)

project(Eigenvalues)

enable_language(CXX)

add_subdirectory("input" ${CMAKE_CURRENT_BINARY_DIR}/inputbin)

add_subdirectory("print" ${CMAKE_CURRENT_BINARY_DIR}/printbin)

add_subdirectory("power_method" ${CMAKE_CURRENT_BINARY_DIR}/power_methodbin)

add_executable(exec main.cpp)

target_include_directories(exec PUBLIC "input")

target_include_directories(exec PUBLIC "print")

target_include_directories(exec PUBLIC "power_method")

target_link_libraries(exec PUBLIC input print power_method)