project(gtests)

FILE(GLOB QUIRC_FILES ../quirc/lib/*)
add_library(quirc ${QUIRC_FILES})
target_include_directories(quirc PUBLIC "../quirc/lib/")

FILE(GLOB QRGEN_FILES ../qrgen/cpp/*)
add_library(qrgen ${QRGEN_FILES})
target_include_directories(quirc PUBLIC "../qrgen/cpp/")

add_executable(
  hello_test
  hello_test.cc
  compression.cpp
  compression.hpp
)

target_link_libraries(
  hello_test
  gtest_main
  qrencode
  quirc
  zbar
  qrgen
)

include(GoogleTest)
gtest_discover_tests(hello_test)
