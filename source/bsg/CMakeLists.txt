set(bsg_sources src/Generator.cc src/OptionContainer.cc src/SpectralFunctions.cc src/Utilities.cc)
set(bsg_headers include/ChargeDistributions.h include/Constants.h include/Generator.h include/OptionContainer.h include/Screening.h include/SpectralFunctions.h include/Utilities.h)

add_library(bsg ${bsg_sources} ${bsg_headers})

target_link_libraries(bsg nme ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_custom_command(TARGET bsg
                   POST_BUILD
                 COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bsg> ${PROJECT_BINARY_DIR}/lib/$<TARGET_FILE_NAME:bsg>)

install(TARGETS bsg EXPORT bsg-targets
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib/bsg)

install(FILES ${bsg_headers} DESTINATION include)
