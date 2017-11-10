IF(DOXYGEN)
    IF(DOXYGEN_DOT_PATH)
        SET(DOXYGEN_HAVE_DOT "YES")
    ELSE(DOXYGEN_DOT_PATH)
        SET(DOXYGEN_HAVE_DOT "NO")
    ENDIF(DOXYGEN_DOT_PATH)

    # Configure the doxygen file
    CONFIGURE_FILE (${SLSM_SOURCE_DIR}/dox/Doxyfile.in ${SLSM_BINARY_DIR}/dox/Doxyfile)

    ADD_CUSTOM_COMMAND (OUTPUT ${SLSM_BINARY_DIR}/doc/html/index.html
                        COMMAND ${DOXYGEN} ${SLSM_BINARY_DIR}/dox/Doxyfile
                        DEPENDS ${SLSM_SOURCE_DIR}/dox/*.dox
                                ${SLSM_SOURCE_DIR}/src/*.h
                                ${SLSM_SOURCE_DIR}/src/*.cpp
                                ${SLSM_SOURCE_DIR}/dox/Doxyfile.in)

    ADD_CUSTOM_TARGET (doc DEPENDS ${SLSM_BINARY_DIR}/doc/html/index.html)
ENDIF (DOXYGEN)
