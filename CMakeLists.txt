cmake_minimum_required(VERSION 3.20)
project(poly C)

set(CMAKE_C_STANDARD 99)

add_executable(poly main.c poly.c poly.h)
