cmake_minimum_required(VERSION 3.10)
project(resistors)

if(CMAKE_GENERATOR MATCHES "Visual Studio")
	include(${CMAKE_BINARY_DIR}/conanbuildinfo_multi.cmake)
else()
	include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
endif()

conan_basic_setup(TARGETS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(resistors 
    formatters.h
    main.cpp 
    math_utils.h
    math_utils.cpp
    model.h
    resistor.h
    resistor_series.h
    resistor_series.cpp
    signals.h
    signals.cpp)

target_link_libraries(resistors PRIVATE CONAN_PKG::fmt)
