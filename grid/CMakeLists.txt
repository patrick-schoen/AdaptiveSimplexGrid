set(GRIDFILES 2dmesh.txt 3dmesh.txt)
foreach(file ${GRIDFILES})
  configure_file(${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach()
