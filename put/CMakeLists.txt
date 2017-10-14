
include_directories(${PROJECT_SOURCE_DIR} ${PPA_Assembler_EXTERNAL_INCLUDES})
link_directories(${PPA_Assembler_EXTERNAL_LIBRARIES})

add_executable(put put.cpp)
target_link_libraries(put PPA_Assembler-iniparser ${COMMON_LINK_LIBS})
