# find_program(PYTHON "python")
#
# if (PYTHON)
#     set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
#     set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/module/__init__.py")
#     set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
#
#     add_custom_command(OUTPUT ${OUTPUT}
#                        COMMAND ${PYTHON} ${SETUP_PY} build
#                        COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
#                        DEPENDS ${DEPS})
#
#     add_custom_target(target ALL DEPENDS ${OUTPUT})
#
#     install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
# endif()
