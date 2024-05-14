include(testing.cmake)

add_executable(test_defragmentation EXCLUDE_FROM_ALL test_defragmentation.cc
  "${CMAKE_SOURCE_DIR}/pcAdapter2.cc"
)

target_include_directories(test_defragmentation PRIVATE "${CMAKE_SOURCE_DIR}")
target_link_libraries(test_defragmentation PRIVATE SCOREC::core)
