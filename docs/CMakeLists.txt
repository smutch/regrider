find_package(Doxygen)
find_package(Sphinx)
if(Doxygen_FOUND AND Sphinx_FOUND)

    get_target_property(REGRIDER_PUBLIC_HEADER_DIR regrider INTERFACE_INCLUDE_DIRECTORIES)
    file(GLOB_RECURSE REGRIDER_PUBLIC_HEADERS ${REGRIDER_PUBLIC_HEADER_DIR}/*.hpp)

    set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
    add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                       DEPENDS ${CAT_CUTIFIER_PUBLIC_HEADERS}
                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                       MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                       COMMENT "Generating docs")
     
    add_custom_target(Doxygen DEPENDS ${DOXYGEN_INDEX_FILE})


     
    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)
     
    add_custom_target(docs DEPENDS Doxygen
                      COMMAND ${SPHINX_EXECUTABLE} -b html
                      # Tell Breathe where to find the Doxygen output
                      -Dbreathe_projects.Regrider=${DOXYGEN_OUTPUT_DIR}/xml
                      ${SPHINX_SOURCE} ${SPHINX_BUILD}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")

else()

    message(WARNING "Failed to find both Doxygen and Sphinx. You will not be able to generate docs.")

endif()
