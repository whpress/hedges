add_library(NRpyDNAcode SHARED NRpyDNAcode.cpp)
target_include_directories(NRpyDNAcode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                          ${PYTHON_27_INCLUDE} 
                                          ${NUMPY_1_13_INCLUDE}
)
add_library(NRpyRS SHARED NRpyRS.cpp)
target_include_directories(NRpyRS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                        ${PYTHON_27_INCLUDE} 
                                        ${NUMPY_1_13_INCLUDE}
                                        ${SRC}
)
set_target_properties(NRpyDNAcode NRpyRS
                      PROPERTIES PREFIX ""
)