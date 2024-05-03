cmake_minimum_required(VERSION 3.27)
project(BIN_projekt)

set(CMAKE_CXX_STANDARD 98)

add_executable(BIN_projekt main.cpp cgp.h symb.h
        cgpROC.h)
add_definitions(-DDONOTEVALUATEUNUSEDNODES)
