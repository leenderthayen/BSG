set(bsg_sources src/Generator.cc src/BSGOptionContainer.cc src/SpectralFunctions.cc src/Utilities.cc)
set(bsg_headers include/ChargeDistributions.h include/Constants.h include/Generator.h include/BSGOptionContainer.h include/Screening.h include/SpectralFunctions.h include/Utilities.h)

add_library(bsg_static STATIC ${bsg_sources})
add_library(bsg SHARED ${bsg_sources})

target_link_libraries(bsg nme ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bsg_static nme_static ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(TARGET bsg
                   POST_BUILD
                 COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bsg> ${PROJECT_BINARY_DIR}/lib/$<TARGET_FILE_NAME:bsg>)

add_custom_command(TARGET bsg_static
                  POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bsg_static> ${PROJECT_BINARY_DIR}/lib/$<TARGET_FILE_NAME:bsg_static>)

install(TARGETS bsg bsg_static EXPORT bsg-targets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

install(FILES ${bsg_headers} DESTINATION include)
