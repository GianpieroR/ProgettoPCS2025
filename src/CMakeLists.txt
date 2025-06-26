list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Triangolazione.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Triangolazione1.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Triangolazione2.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Polyhedra.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Dimensione.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/MeshExport.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Duale.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/CamminoMinimo.cpp)

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_headers ${src_headers} PARENT_SCOPE)
set(src_sources ${src_sources} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)