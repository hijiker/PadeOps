file(GLOB boltproblems_source_files *.F90)

# Include directories
include_directories( ${FFTW_INCLUDE_PATH} ${DECOMP_INCLUDE_PATH} ${PadeOps_BINARY_DIR}/src  )

# Link directories
link_directories( ${FFTW_LIBRARY_PATH} ${DECOMP_LIBRARY_PATH} ${PadeOps_BINARY_DIR}/src )

# Create the bolt problem executables
foreach ( boltproblemfile ${boltproblems_source_files} )
    string( REPLACE "${PadeOps_SOURCE_DIR}/problems/bolt/" "" boltproblemfile ${boltproblemfile} )
    string( REPLACE ".F90" "" boltproblemexec ${boltproblemfile} )
  
    add_executable( ${boltproblemexec} ${boltproblemfile} )
    
    target_link_libraries( ${boltproblemexec} BoltOps 2decomp_fft fftw3)

    if (MPI_Fortran_COMPILER_FLAGS)
        set_target_properties(${boltproblemexec} PROPERTIES COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
    
    file(MAKE_DIRECTORY ${PadeOps_BINARY_DIR}/problems/bolt/${boltproblemexec}_files)
    file(COPY ${boltproblemexec}_files/input.dat DESTINATION ${PadeOps_BINARY_DIR}/problems/bolt/${boltproblemexec}_files) 
endforeach()
