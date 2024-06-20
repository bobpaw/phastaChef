include(testing.cmake)

add_executable(test_defragmentation EXCLUDE_FROM_ALL test_defragmentation.cc
  "${CMAKE_SOURCE_DIR}/pcAdapter2.cc"
)

target_include_directories(test_defragmentation PRIVATE "${CMAKE_SOURCE_DIR}")
target_link_libraries(test_defragmentation PRIVATE SCOREC::core)

add_executable(test_raytracing EXCLUDE_FROM_ALL test_raytracing.cc
  "${CMAKE_SOURCE_DIR}/pcAdapter2.cc"
)
target_link_libraries(test_raytracing PRIVATE SCOREC::core)

add_executable(test_labeling EXCLUDE_FROM_ALL test_labeling.cc
  "${CMAKE_SOURCE_DIR}/pcAdapter2.cc"
)

target_include_directories(test_labeling PRIVATE "${CMAKE_SOURCE_DIR}")
target_link_libraries(test_labeling PRIVATE SCOREC::core)

add_executable(test_prisms EXCLUDE_FROM_ALL test_prisms.cc
  "${CMAKE_SOURCE_DIR}/pcAdapter2.cc"
)

target_include_directories(test_prisms PRIVATE "${CMAKE_SOURCE_DIR}")
target_link_libraries(test_prisms PRIVATE SCOREC::core)
