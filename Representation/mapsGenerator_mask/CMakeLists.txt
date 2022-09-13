project(maps_generator)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

if(UNIX AND NOT APPLE)
list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -Wall -pthread -DSCWRL4_ENERGY -O3")
else()
list(APPEND CMAKE_CXX_FLAGS "-std=c++0x -Wall -DSCWRL4_ENERGY -O3")
endif()

set(CMAKE_VERBOSE_MAKEFILE 1)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)

cmake_policy(SET CMP0015 NEW)
link_directories(
)

cmake_minimum_required(VERSION 2.6.2)

include_directories(
  libs
  libs/proteins
  libs/c2numpy
  libs/libnpy
  /usr/local/include/eigen3/
)

add_subdirectory(
  libs/proteins
)

add_executable(maps_generator
  src/engine.cpp
  src/main.cpp
)

if(UNIX AND NOT APPLE)
    #SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -static -lc++" )
    SET( cmake_exe_linker_flags  "${cmake_exe_linker_flags} -static -lc++" )
    target_link_libraries( maps_generator proteins  )
elseif(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "-static  -lc++")
  target_link_libraries( maps_generator proteins)
else()

    SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -static -lc++" )
    SET( cmake_exe_linker_flags  "${cmake_exe_linker_flags} -static -lc++" )
    target_link_libraries( maps_generator proteins  )

    # Mac OS X specific post code
    add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    COMMAND strip Release/${PROJECT_NAME}
)

endif()

