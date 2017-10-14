
include_directories(${PROJECT_SOURCE_DIR} ${PPA_Assembler_EXTERNAL_INCLUDES})
link_directories(${PPA_Assembler_EXTERNAL_LIBRARIES})

add_library(iniparser-objs OBJECT
    dictionary.c
    iniparser.c
    )

add_library(PPA_Assembler-iniparser STATIC $<TARGET_OBJECTS:iniparser-objs>)

target_link_libraries(PPA_Assembler-iniparser ${COMMON_LINK_LIBS})
target_include_directories (PPA_Assembler-iniparser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
