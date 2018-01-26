add_executable(bsg_exec BSG.cc)

target_link_libraries(bsg_exec bsg)

install(TARGETS bsg_exec EXPORT BSGTargets
	RUNTIME DESTINATION bin)
