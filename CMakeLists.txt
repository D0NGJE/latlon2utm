cmake_minimum_required(VERSION 3.0.2)
project(latlon_to_utm)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  nav_msgs
  roscpp
  sensor_msgs
  std_msgs
)

include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)

add_executable(latlon2utm src/latlon2utm.cpp)
target_link_libraries(latlon2utm
  ${catkin_LIBRARIES}
)
