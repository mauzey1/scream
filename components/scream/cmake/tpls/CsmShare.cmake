macro (CreateCsmShareTarget)

  # Some sanity checks
  if (NOT SCREAM_CIME_BUILD)
    message (FATAL_ERROR "Error! CsmShare.cmake currently only works in a CIME build.")
  endif ()
  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()
  if (NOT DEFINED COMP_INTERFACE)
    message (FATAL_ERROR "Error! The cmake variable 'COMP_INTERFACE' is not defined.")
  endif ()
  if (NOT DEFINED NINST_VALUE)
    message (FATAL_ERROR "Error! The cmake variable 'NINST_VALUE' is not defined.")
  endif ()

  if (USE_ESMF_LIB)
    set(ESMFDIR "esmf")
  else()
    set(ESMFDIR "noesmf")
  endif()

  set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")

  # If we didn't already parse this script, create interface lib
  if (NOT TARGET scream_csm_share)
    find_library(CSM_SHARE_LIB csm_share REQUIRED PATHS ${CSM_SHARE})

    # Create the imported library that scream targets can link to
    add_library(scream_csm_share UNKNOWN IMPORTED GLOBAL)
    set_target_properties(scream_csm_share PROPERTIES IMPORTED_LOCATION ${CSM_SHARE_LIB})
    set_target_properties(scream_csm_share PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CSM_SHARE})

    include(tpls/Scorpio)
    CreateScorpioTarget(TRUE)
    target_link_libraries(scream_csm_share INTERFACE scream_piof)
  endif ()
endmacro()
