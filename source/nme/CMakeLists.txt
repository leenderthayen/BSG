set(nme_sources src/NMEOptions.cc src/NuclearStructureManager.cc)
set(nme_headers include/MatrixElements.h include/NilssonOrbits.h include/NMEOptions.h include/NuclearStructureManager.h include/NuclearUtilities.h)

add_library(nme ${nme_sources} ${nme_headers})

target_link_libraries(nme ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(TARGET nme
                   POST_BUILD
                 COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:nme> ${PROJECT_BINARY_DIR}/lib/$<TARGET_FILE_NAME:nme>)

install(TARGETS nme EXPORT bsg-targets
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib/nme)

install(FILES ${nme_headers} DESTINATION include)
