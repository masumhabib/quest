
# Automatic versioning using git and static file.

# 
# Returns the version string from git if available. Ohterwise, get it from
# the VERSION file.
#
set(VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
function(get_version major minor patch build)
    git_describe (version_from_git)

    if (EXISTS "${VERSION_FILE}")
        file (STRINGS "${VERSION_FILE}" version_from_file)
    else ()
        set (version_from_file "v0.0.0")
    endif()

    if (version_from_git MATCHES ".*-NOTFOUND")
        set(version ${version_from_file})
    else ()
        set(version ${version_from_git})
    endif()
    
    #message (STATUS "DBG: ${version}")

    string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" _major "${version}")
    string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" _minor "${version}")
    string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" _patch "${version}")
    string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" _build "${version}")

    # update the VERSION file if git version does not match with it.
    if (NOT "v${_major}.${_minor}.${_patch}" MATCHES "^${version_from_file}$")
        file(WRITE ${VERSION_FILE} "v${_major}.${_minor}.${_patch}")
    endif()


    set(${major} ${_major} PARENT_SCOPE) 
    set(${minor} ${_minor} PARENT_SCOPE) 
    set(${patch} ${_patch} PARENT_SCOPE) 
    set(${build} ${_build} PARENT_SCOPE) 
endfunction()

# 
# Returns the version string from Git tags
#
# Inspired from 
# https://github.com/rpavlik/cmake-modules/blob/master/GetGitRevisionDescription.cmake
#
function(git_describe description)
	if(NOT GIT_FOUND)
		find_package(Git QUIET)
	endif()
	if(NOT GIT_FOUND)
		set(${description} "GIT-NOTFOUND" PARENT_SCOPE)
		return()
	endif()

	execute_process(COMMAND
		"${GIT_EXECUTABLE}" describe
		WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE status
        OUTPUT_VARIABLE out
		ERROR_QUIET 
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

	if(NOT status EQUAL 0)
		set(out "${out}-${status}-NOTFOUND")
	endif()

	set(${description} "${out}" PARENT_SCOPE)
endfunction()




