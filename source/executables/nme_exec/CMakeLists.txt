add_executable(nme_exec NME.cc)

target_link_libraries(nme_exec nme)

add_custom_command(TARGET nme_exec
                   POST_BUILD
                 COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:nme_exec> ${PROJECT_BINARY_DIR}/bin/$<TARGET_FILE_NAME:nme_exec>)

install(TARGETS nme_exec EXPORT bsg-targets
        RUNTIME DESTINATION bin)
