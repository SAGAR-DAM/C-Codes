cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Set the project name
project(Libtorch)

# Set the installation prefix (modify this as needed)
set(CMAKE_INSTALL_PREFIX "D:/C++modules/libtorch/install")

# Find the Torch package
find_package(Torch REQUIRED)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# Add a check to see if we're using MinGW
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -Wa,-mbig-obj")
elseif (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /bigobj")
  file(GLOB TORCH_DLLS "${TORCH_INSTALL_PREFIX}/lib/*.dll")
  add_custom_command(TARGET Libtorch
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${TORCH_DLLS}
                     $<TARGET_FILE_DIR:Libtorch>)
endif()

# Make sure MSVC-specific flags are removed
string(REPLACE "/EHsc" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REPLACE "/bigobj" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Define the executable
add_executable(Libtorch main.cpp)

# Link the Torch library
target_link_libraries(Libtorch "${TORCH_LIBRARIES}")

# Set C++ standard
set_property(TARGET Libtorch PROPERTY CXX_STANDARD 17)

