# Project documentation builder

if(DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/doc/Doxyfile ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/Doxyfile @ONLY)

    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/Doxyfile
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

    add_custom_command(
        TARGET doc  POST_BUILD
        COMMAND mkdir -p ${QUEST_DIST_DIR}/doc
        COMMAND cp -rf ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/html ${QUEST_DIST_DIR}/doc 
        COMMAND cd ${QUEST_DIST_DIR}/doc && ln -sf html/index.html index.html
    )
endif(DOXYGEN_FOUND)


