# message("SUBFOLDER:")

message("CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
message("CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")
message("CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")
message("CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
message("PROJECT_NAME = ${PROJECT_NAME}")
message("CMAKE_PROJECT_NAME = ${CMAKE_PROJECT_NAME}")

list(APPEND src_test_headers "${CMAKE_CURRENT_SOURCE_DIR}/test.hpp")

list(APPEND src_includes "${CMAKE_CURRENT_SOURCE_DIR}")

set(src_test_headers ${src_headers} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)