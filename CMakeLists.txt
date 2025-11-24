cmake_minimum_required(VERSION 3.16)

project(Gamma3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Gamma3 main.cpp)

target_include_directories(Gamma3  PRIVATE
        ./include
        ./third_party

)
