IF(DOXYGEN)
    IF(DOXYGEN_DOT_PATH)
        SET(DOXYGEN_HAVE_DOT "YES")
    ELSE(DOXYGEN_DOT_PATH)
        SET(DOXYGEN_HAVE_DOT "NO")
    ENDIF(DOXYGEN_DOT_PATH)

    # Configure the doxygen file
    CONFIGURE_FILE (${LSM_SOURCE_DIR}/dox/Doxyfile.in ${LSM_BINARY_DIR}/dox/Doxyfile)

    ADD_CUSTOM_COMMAND (OUTPUT ${LSM_BINARY_DIR}/doc/html/index.html
                        COMMAND ${DOXYGEN} ${LSM_BINARY_DIR}/dox/Doxyfile
                        DEPENDS #${LSM_SOURCE_DIR}/dox/*.dox
                                ${LSM_SOURCE_DIR}/src/*.h
                                #${LSM_SOURCE_DIR}/src/*.cpp
                                ${LSM_SOURCE_DIR}/dox/Doxyfile.in)

    ADD_CUSTOM_TARGET (doc DEPENDS ${LSM_BINARY_DIR}/doc/html/index.html)
ENDIF (DOXYGEN)
