macro (CreateMctTarget)

  # Some sanity checks
  if (NOT SCREAM_CIME_BUILD)
    message (FATAL_ERROR "Error! Mct.cmake currently only works in a CIME build.")
  endif ()
  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If we didn't already parse this script, create interface lib
  if (NOT TARGET scream_mct)
    message ("looking for mct lib in ${INSTALL_SHAREDPATH}/lib")
    find_library(MCT_LIB mct REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the imported library that scream targets can link to
    add_library(scream_mct UNKNOWN IMPORTED GLOBAL)
    set_target_properties(scream_mct PROPERTIES IMPORTED_LOCATION ${MCT_LIB})
    set_target_properties(scream_mct PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)

    include(cime/CsmShare)
    CreateCsmShareTarget()
    target_link_libraries(scream_mct INTERFACE scream_csm_share)
  endif ()
endmacro()
