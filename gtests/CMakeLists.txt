project(gtests)

add_executable(
  hello_test
  hello_test.cc
)

target_link_libraries(
  hello_test
  gtest_main
  qrencode
)

include(GoogleTest)
gtest_discover_tests(hello_test)
