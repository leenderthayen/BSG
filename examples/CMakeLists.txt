project(BSG-examples CXX)

find_package(GSL REQUIRED)

set(BSG_EXAMPLES
standalone/bsg_standalone.cc)

foreach(EXAMPLE_NAME IN ITEMS ${BSG_EXAMPLES})
  get_filename_component(EXAMPLE_NAME_NOEXT ${EXAMPLE_NAME} NAME_WE)
  add_executable(${EXAMPLE_NAME_NOEXT} ${EXAMPLE_NAME})
  target_include_directories(${EXAMPLE_NAME_NOEXT} PUBLIC ${GSL_INCLUDE_DIRS})
  target_link_libraries(${EXAMPLE_NAME_NOEXT} PRIVATE BSG)
  target_link_libraries(${EXAMPLE_NAME_NOEXT} PUBLIC ${GSL_LIBRARIES})
endforeach()

# add_custom_command(TARGET bsg_exec
#                    POST_BUILD
#                  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bsg_exec> ${PROJECT_BINARY_DIR}/bin/$<TARGET_FILE_NAME:bsg_exec>)
#
# add_custom_command(TARGET bsg_exec POST_BUILD
#                   COMMAND ${CMAKE_COMMAND} -E copy_directory
#                 ${CMAKE_SOURCE_DIR}/data ${PROJECT_BINARY_DIR}/bin/)
#
# install(TARGETS bsg_exec EXPORT bsg-targets
# 	RUNTIME DESTINATION bin)
