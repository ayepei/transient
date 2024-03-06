add_subdirectory(moc)


add_executable(OPENMOC-transient OPENMOC-transient.cpp )
target_link_libraries(OPENMOC-transient PUBLIC ${MOC_LIBS} )
target_include_directories(OPENMOC-transient PUBLIC ${PROJECT_BINARY_DIR} ${MOC_INCLUDE_DIR} )

add_executable(2OPENMOC-transient OPENMOC-transient2.cpp )
target_link_libraries(2OPENMOC-transient PUBLIC ${MOC_LIBS} )
target_include_directories(2OPENMOC-transient PUBLIC ${PROJECT_BINARY_DIR} ${MOC_INCLUDE_DIR} )
