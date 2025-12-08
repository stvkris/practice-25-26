cmake_minimum_required(VERSION 3.10)

project(LSSolution)

enable_language(CXX)



find_package(Threads REQUIRED)

add_subdirectory("input" ${CMAKE_CURRENT_BINARY_DIR}/inputbin)

add_subdirectory("print" ${CMAKE_CURRENT_BINARY_DIR}/printbin)

add_subdirectory("solution" ${CMAKE_CURRENT_BINARY_DIR}/solutionbin)

add_executable(exec main.cpp)

target_include_directories(exec PUBLIC "input")

target_include_directories(exec PUBLIC "print")

target_include_directories(exec PUBLIC "solution")

target_link_libraries(exec PUBLIC input print solution)

target_link_libraries(exec PUBLIC Threads::Threads)

