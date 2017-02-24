# create a plugin that adds a reader to the ParaView GUI# it is added in the file dialog when doing opens/saves.

# The qrc file is processed by Qt's resource compiler (rcc)
# the qrc file must have a resource prefix of "/ParaViewResources"
# and ParaView will read anything contained under that prefix
# the pqReader.xml file contains xml defining readers with their
# file extensions and descriptions.

INCLUDE(${PARAVIEW_USE_FILE})

# Set output library destination to plugin folder
#SET(
#  LIBRARY_OUTPUT_PATH $ENV{PV_PLUGIN_PATH}
#  CACHE INTERNAL
#  "Single output directory for building all libraries."
#)

# Build the server-side plugin
ADD_PARAVIEW_PLUGIN(
  PV4FoamReader_SM "1.0"
  SERVER_MANAGER_XML PV4FoamReader_SM.xml
  SERVER_MANAGER_SOURCES vtkPV4FoamReader.cxx
)

TARGET_LINK_LIBRARIES(PV4FoamReader_SM PUBLIC vtkPV4Foam pqCore)

#-----------------------------------------------------------------------------
