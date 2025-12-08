cmake_minimum_required(VERSION 3.10)

project(LSSSolution)

enable_language(CXX)

find_package(Threads REQUIRED)

add_library(solution STATIC solution.cpp)

target_link_libraries(solution PUBLIC Threads::Threads)
