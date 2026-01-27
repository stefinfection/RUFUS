SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(JELLYFISH_PROJECT jellyfish_project CACHE INTERNAL "jellyfish project name")
SET(JELLYFISH_DIR ${CMAKE_BINARY_DIR}/externals/jellyfish CACHE INTERNAL "jellyfish project directory")
SET(JELLYFISH_LIB)

ExternalProject_Add(${JELLYFISH_PROJECT}
	GIT_REPOSITORY https://github.com/gmarcais/Jellyfish.git
	GIT_TAG master
	CONFIGURE_COMMAND "pwd; ./configure --prefix=/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/bin/externals/jellyfish/src/jellyfish_project-build/"
	BUILD_IN_SOURCE 1
	#CONFIGURE_COMMAND ./configure
	BUILD_COMMAND make
	INSTALL_COMMAND make install
	UPDATE_COMMAND ""
	PREFIX ${JELLYFISH_DIR}
)

# For building on newer system
# ExternalProject_Add(${MODIFIED_JELLYFISH_PROJECT}
# 	URL ${PROJECT_SOURCE_DIR}/src/modifiedJellyfish.tar.gz

#  	PATCH_COMMAND
#     		patch -p0 < ${PROJECT_SOURCE_DIR}/src/externals/patches/jellyfish-configure-ac.patch

#  	BUILD_IN_SOURCE 1

# 	CONFIGURE_COMMAND
#     		autoreconf -fi
#     		COMMAND ./configure --prefix=${PROJECT_SOURCE_DIR}/bin/externals/modified_jellyfish/src/modified_jellyfish_project/

#   	BUILD_COMMAND make
#         INSTALL_COMMAND make install
# 	#UPDATE_COMMAND ""
# 	#PREFIX ${MODIFIED_JELLYFISH_DIR}
# )

ExternalProject_Get_Property(${JELLYFISH_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${JELLYFISH_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${JELLYFISH_PROJECT} BINARY_DIR)

SET(JELLYFISH_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "JELLYFISH INCLUDE")