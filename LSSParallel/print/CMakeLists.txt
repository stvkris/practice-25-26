cmake_minimum_required(VERSION 3.10)

project(LSSPrint)

enable_language(CXX)

find_package(Threads REQUIRED)

add_library(print STATIC print.cpp)

target_link_libraries(print PUBLIC Threads::Threads)
