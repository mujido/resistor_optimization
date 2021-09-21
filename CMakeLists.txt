cmake_minimum_required(VERSION 3.10)

project(resistors)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(resistors 
    main.cpp 
    math_utils.h
    math_utils.cpp
    resistor.h
    resistor_series.h
    resistor_series.cpp
    signals.h
    signals.cpp)

if(MSVC)
   target_include_directories(resistors PUBLIC 
        "C:\\Users\\kurt\\Documents\\C++\\ThirdParty\\boost_1_70_0"
    )
endif()