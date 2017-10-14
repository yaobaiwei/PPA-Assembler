
include_directories(${PROJECT_SOURCE_DIR} ${PPA_Assembler_EXTERNAL_INCLUDES})
link_directories(${PPA_Assembler_EXTERNAL_LIBRARIES})

#### Main ####
add_executable(run example.cpp)
target_link_libraries(run PPA_Assembler-iniparser ${COMMON_LINK_LIBS})
