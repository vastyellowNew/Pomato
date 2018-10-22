set(MAPS_SRC
	newton.cpp
)
add_library(maps ${LIBMODE}
	${MAPS_SRC}
)
target_link_libraries( maps ${Teem_LIBRARIES} )

install( TARGETS maps
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

set(MAPS_EXEC_SRC
	map_ftle
	parallel_analysis_stdmap
	test_stdmap
)

foreach( target ${MAPS_EXEC_SRC} )
  	add_executable( ${target} ${target}.cpp )
	target_link_libraries( ${target} maps 
        ${GLUT_LIBRARIES} 
        ${OPENGL_LIBRARIES}
        ${GLUI_LIBRARIES}
        ${Teem_LIBRARIES} )
	install( TARGETS ${target} RUNTIME DESTINATION bin )
endforeach( target ${MAPS_EXEC_SRC} )
