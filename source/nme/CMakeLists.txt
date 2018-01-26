set(nme_sources src/NMEOptions.cc src/NuclearStructureManager.cc)
set(nme_headers include/MatrixElements.h include/NilssonOrbits.h include/NMEOptions.h include/NuclearStructureManager.h include/NuclearUtilities.h)

add_library(nme ${nme_sources})

target_link_libraries(nme ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS nme EXPORT BSGTargets
        RUNTIME DESTINATION bin 
        LIBRARY DESTINATION lib 
        ARCHIVE DESTINATION lib/nme)

install(FILES ${nme_headers} DESTINATION include)
