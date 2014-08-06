cmake_minimum_required(VERSION 2.8.3)
project(force_estimation)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  cisst_msgs
  cisst_ros_bridge
  roscpp
  std_msgs
)

set (REQUIRED_CISST_LIBRARIES
     cisstCommon
     cisstVector
     cisstOSAbstraction
     cisstMultiTask
     cisstParameterTypes
     cisstRobot
     cisstNumerical 
     )

# find cisst and make sure required libraries have been compiled
find_package (cisst REQUIRED ${REQUIRED_CISST_LIBRARIES})


catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES force_estimation
#  CATKIN_DEPENDS cisst_msgs cisst_ros_bridge roscpp std_msgs
#  DEPENDS system_lib
)

 
include_directories(
  ${catkin_INCLUDE_DIRS}
)
 if( cisst_FOUND)

 include (${CISST_USE_FILE})

  cisst_find_saw_component( sawCANBus REQUIRED )
  cisst_find_saw_component( sawBarrett REQUIRED )
  cisst_find_saw_component( sawControllers REQUIRED )
  cisst_find_saw_component( sawKeyboard REQUIRED )
  cisst_find_saw_component( sawJR3ForceSensor REQUIRED )
 
include_directories(
    ${sawCANBus_INCLUDE_DIR}
    ${sawBarrett_INCLUDE_DIR}
    ${sawKeyboard_INCLUDE_DIR}
    ${sawControllers_INCLUDE_DIR}
    ${sawJR3ForceSensor_INCLUDE_DIR} )

  
  link_directories( 
    ${sawCANBus_LIBRARY_DIR}
    ${sawBarrett_LIBRARY_DIR}
    ${sawControllers_LIBRARY_DIR} 
    ${sawKeyboard_LIBRARY_DIR}
    ${sawJR3ForceSensor_LIBRARY_DIR} )

  
  set(
    REQUIRED_SAW_COMPONENTS
    sawCANBus
    sawBarrett
    sawControllers
    sawJR3ForceSensor
    sawKeyboard )

add_executable( 
      wamForceEstimation 
      src/wamForceEstimation.cpp 
      src/osaHybridForcePosition.cpp
      src/mtsHybridForcePosition.cpp )

  target_link_libraries( wamForceEstimation ${REQUIRED_SAW_COMPONENTS} ${catkin_LIBRARIES})
  cisst_target_link_libraries( wamForceEstimation ${REQUIRED_CISST_LIBRARIES})
  
endif( cisst_FOUND )
## Generate messages in the 'msg' folder
# add_message_files(
#   FILES
#   Message1.msg
#   Message2.msg
# )

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
# generate_messages(
#   DEPENDENCIES
#   cisst_msgs#   std_msgs
# )



