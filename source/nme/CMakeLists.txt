set(nme_sources src/NMEOptionContainer.cc src/NuclearStructureManager.cc)
set(nme_headers include/MatrixElements.h include/NilssonOrbits.h include/NMEOptionContainer.h include/NuclearStructureManager.h include/NuclearUtilities.h)

add_library(nme_static STATIC ${nme_sources})
add_library(nme SHARED ${nme_sources})

target_link_libraries(nme ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(nme_static ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(TARGET nme
                   POST_BUILD
                 COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:nme> ${PROJECT_BINARY_DIR}/lib/$<TARGET_FILE_NAME:nme>)
add_custom_command(TARGET nme_static
                  POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:nme_static> ${PROJECT_BINARY_DIR}/lib/$<TARGET_FILE_NAME:nme_static>)

install(TARGETS nme nme_static EXPORT bsg-targets
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib/nme
ARCHIVE DESTINATION lib/nme)

install(FILES ${nme_headers} DESTINATION include)
