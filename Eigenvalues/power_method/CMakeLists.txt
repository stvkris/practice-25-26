cmake_minimum_required(VERSION 3.10)

project(EigenPower)

enable_language(CXX)

add_library(power_method STATIC power_method.cpp)
