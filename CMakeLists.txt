# $Id: CMakeLists.txt Majid Kazemi$
# Purpose. Event Identification $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(EI_Project)

#SET(CMAKE_INSTALL_PREFIX "/home/majid/SiFi-CC/EventIdenfication/build1")
SET(CMAKE_INSTALL_PREFIX "/home/majid/EveIdenMVA/build")
message(STATUS "Install path " ${CMAKE_INSTALL_PREFIX})

#----------Set Compiler
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS} ${CMAKE_CURRENT_SOURCE_DIR} )
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/majid/root-git/root-build" ${CMAKE_CURRENT_SOURCE_DIR} )
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/majid/root6.13/root-build" ${CMAKE_CURRENT_SOURCE_DIR} )

message(STATUS "Found CXXCompiler version: " ${CMAKE_CXX_COMPILER_VERSION})

#--------------Root Setup
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#list(APPEND CMAKE_PREFIX_PATH "/home/majid/root-git/root-build")
#list(APPEND CMAKE_PREFIX_PATH "/home/majid/root6.13/root-build")

find_package( ROOT REQUIRED COMPONENTS Gui Core) 
message( STATUS "Found ROOT version " ${ROOT_VERSION} ) 

include(${ROOT_USE_FILE})

include(makeRootDict)


include_directories(${ROOT_INCLUDE_DIRS} )

# Setup include directory for this project

include_directories(${PROJECT_SOURCE_DIR}/include)

#---------------------------------------------------------------------------
#create extra libraries

MAKE_ROOT_DICT(DICTNAME StructLib LINKDEF include/LinkDef.h SOURCES src/DR_GenerallStructs.cc INCLUDES include/Structinclude.hh include/DR_GenerallStructs.hh LIBRARIES ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries


SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib )


#add_executable(EventIdenfication EventIdenfication.cc ${sources} ${headers})
#target_link_libraries(EventIdenfication StructLib ${ROOT_LIBRARIES})

#add_executable(EveIdenMVA EveIdenMVA.cc ${sources} ${headers})
#target_link_libraries(EveIdenMVA StructLib ${ROOT_LIBRARIES})

add_executable(CompId CompId.cc ${sources} ${headers})
target_link_libraries(CompId StructLib ${ROOT_LIBRARIES})

#add_executable(ICPS ICPS.cc ${sources} ${headers})
#target_link_libraries(ICPS StructLib ${ROOT_LIBRARIES})

#---Create  a main program using the library

#target_link_libraries(EventIdenfication EI_Project)


#----------------------------------------------------------------------------

#add_custom_target(EI_Project DEPENDS EventIdenfication)
#add_custom_target(EI_Project DEPENDS EveIdenMVA)
add_custom_target(EI_Project DEPENDS CompId)
#add_custom_target(EI_Project DEPENDS ICPS)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX

#set_property(TARGET EventIdenfication PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#install(TARGETS EventIdenfication RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})

#set_property(TARGET EveIdenMVA PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#install(TARGETS EveIdenMVA RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})

set_property(TARGET CompId PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS CompId RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})

#set_property(TARGET ICPS PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#install(TARGETS ICPS RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
