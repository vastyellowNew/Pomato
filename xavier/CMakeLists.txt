# individual directorie
list( APPEND XAVIER_DIRECTORIES cr3bp )
list( APPEND XAVIER_DIRECTORIES convolution )
list( APPEND XAVIER_DIRECTORIES vis )

foreach( lib ${XAVIER_DIRECTORIES})
  add_subdirectory( ${lib} )
endforeach()

######################################################################
# determine targets to build and define their linking
######################################################################

install( DIRECTORY bin DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY scripts DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY utils/GLSLShaderCode DESTINATION utils USE_SOURCE_PERMISSIONS)

######################################################################
# add local libraries to targets
######################################################################
# foreach( target ${XAVIER_TARGETS} )
#     list( APPEND PMATE_TARGETS $target )
# endforeach()


