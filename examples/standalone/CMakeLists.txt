add_executable(bsg_exec BSG.cc)

target_link_libraries(bsg_exec bsg ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(TARGET bsg_exec
                   POST_BUILD
                 COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bsg_exec> ${PROJECT_BINARY_DIR}/bin/$<TARGET_FILE_NAME:bsg_exec>)

add_custom_command(TARGET bsg_exec POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/data ${PROJECT_BINARY_DIR}/bin/)

install(TARGETS bsg_exec EXPORT bsg-targets
	RUNTIME DESTINATION bin)
