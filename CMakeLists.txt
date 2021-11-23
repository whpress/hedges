project(native C CXX)
cmake_minimum_required(VERSION 3.5)
get_filename_component(NUMPY_1_13_INCLUDE /usr/local/lib/python2.7/dist-packages/numpy-1.13.3-py2.7-linux-x86_64.egg/numpy/core/include ABSOLUTE)
get_filename_component(PYTHON_27_INCLUDE /usr/local/include/python2.7/ ABSOLUTE)
get_filename_component(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src ABSOLUTE)
get_filename_component(SCHIFRA ${CMAKE_CURRENT_SOURCE_DIR}/src/schifra ABSOLUTE)

#get_filename_component(PYTHONSRC ${CMAKE_CURRENT_SOURCE_DIR}/python ABSOLUTE)
#get_filename_component(CPYTHONSRC ${CMAKE_CURRENT_SOURCE_DIR}/cpython ABSOLUTE)

add_subdirectory(src)

enable_testing(TRUE)