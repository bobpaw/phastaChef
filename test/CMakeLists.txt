include(testing.cmake)

add_executable(test_labeling test_labeling.cc
  "${CMAKE_SOURCE_DIR}/pcAdapter2.cc"
)

target_include_directories(test_labeling PRIVATE "${CMAKE_SOURCE_DIR}")
target_link_libraries(test_labeling PRIVATE SCOREC::core)
