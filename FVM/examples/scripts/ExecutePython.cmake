cmake_minimum_required(VERSION 3.21)


function(make_python_import_path IMPORT_DIRS)
    if(WIN32)
        set(PATH_SEPARATOR ";")
    else()
        set(PATH_SEPARATOR ":")
    endif()

    foreach(DIR ${IMPORT_DIRS})
        if (NOT $ENV{PYTHONPATH})
            set(ENV{PYTHONPATH} "${DIR}")
        else()
            set(ENV{PYTHONPATH} "$ENV{PYTHONPATH}${PATH_SEPARATOR}${DIR}")
        endif()
    endforeach()
    
endfunction()


function(execute_python_script)
    set(options "")
    set(oneValueArgs SCRIPT_PATH)
    set(multiValueArgs IMPORT_DIRS)
    cmake_parse_arguments("EPS" "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT EPS_SCRIPT_PATH)
        message(FATAL_ERROR "SCRIPT_PATH is required")
    endif()

    set(ORIGINAL_PYTHON_PATH $ENV{PYTHONPATH})
    make_python_import_path(${EPS_IMPORT_DIRS})

    execute_process(
        COMMAND ${Python_EXECUTABLE} 
        ${EPS_SCRIPT_PATH}
        RESULT_VARIABLE PYTHON_RESULT
        OUTPUT_VARIABLE PYTHON_STDOUT
        ERROR_VARIABLE  PYTHON_STDERR
    )

    set(ENV{PYTHONPATH} ${ORIGINAL_PYTHON_PATH})

    if (NOT ${PYTHON_RESULT} EQUAL 0)
        message(FATAL_ERROR
            "Python script returned non zero status ${PYTHON_RESULT}\n"
            ${PYTHON_STDERR}
        )
    endif()

endfunction()
