project(BSG-utests CXX)

find_package(GSL REQUIRED)

set(BSG_MAIN main.cc)

set(BSG_UTESTS_SOURCES
  test_gamma_functions.cc
  test_spectral_functions.cc
  test_config_parser.cc)

enable_testing()

foreach(TEST_NAME IN ITEMS ${BSG_UTESTS_SOURCES})
  get_filename_component(TEST_NAME_NOEXT ${TEST_NAME} NAME_WE)
  add_executable(${TEST_NAME_NOEXT} ${BSG_MAIN} ${TEST_NAME})
  target_include_directories(${TEST_NAME_NOEXT} PUBLIC ${GSL_INCLUDE_DIRS})
  target_link_libraries(${TEST_NAME_NOEXT} PRIVATE BSG)
  target_link_libraries(${TEST_NAME_NOEXT} PUBLIC ${GSL_LIBRARIES})
  add_test(NAME ${TEST_NAME_NOEXT} COMMAND ${TEST_NAME_NOEXT} -r console)
endforeach()
