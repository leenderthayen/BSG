find_program(PYTHON "python")

if (PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/bsg_gui/__init__.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${PYTHON} -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt
                       COMMAND ${PYTHON} ${SETUP_PY} build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/bsg_gui ${PROJECT_BINARY_DIR}/bin/bsg_gui/
                       DEPENDS ${DEPS})

    add_custom_target(bsg_gui ALL DEPENDS ${OUTPUT})

    # Uncomment if you want to actually install it as a python package
    #install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")

    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bsg_gui DESTINATION bin)
endif()
