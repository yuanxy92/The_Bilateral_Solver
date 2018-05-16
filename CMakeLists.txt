# CMake file for GigaRender project
# Author: Shane Yuan
# E-mail: yuanxy92@gmail.com
#

##
##                     _oo0oo_
##                    o8888888o
##                    88" . "88
##                    (| -_- |)
##                    0\  =  /0
##                  ___/`---'\___
##                .' \\|     |// '.
##               / \\|||  :  |||// \
##              / _||||| -:- |||||- \
##             |   | \\\  -  /// |   |
##             | \_|  ''\---/''  |_/ |
##             \  .-\__  '-'  ___/-. /
##           ___'. .'  /--.--\  `. .'___
##        ."" '<  `.___\_<|>_/___.' >' "".
##       | | :  `- \`.;`\ _ /`;.`/ - ` : | |
##       \  \ `_.   \_ __\ /__ _/   .-` /  /
##   =====`-.____`.___ \_____/___.-`___.-'=====
##                     `=---='
##
##             佛祖保佑         永无BUG
##               Buddha Bless No Bug
##

message( 
"##
##                     _oo0oo_
##                    o8888888o
##                    88\" . \"88
##                    (| -_- |)
##                    0\\  =  /0
##                  ___/`---'\\___
##                .' \\\\|     |// '.
##               / \\\\|||  :  |||// \\
##              / _||||| -:- |||||- \\
##             |   | \\\\\\  -  /// |   |
##             | \\_|  ''\\---/''  |_/ |
##             \\  .-\\__  '-'  ___/-. /
##           ___'. .'  /--.--\\  `. .'___
##        .\"\" '<  `.___\\_<|>_/___.' >' \"\".
##       | | :  `- \\`.;`\\ _ /`;.`/ - ` : | |
##       \\  \\ `_.   \\_ __\\ /__ _/   .-` /  /
##   =====`-.____`.___ \\_____/___.-`___.-'=====
##                     `=---='
##
##             佛祖保佑         永无BUG
##               Buddha Bless No Bug
##
"
)

cmake_minimum_required (VERSION 3.10)
set(name "FastBilateralSolver")
project ( ${name} )

# set output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/bin/Debug)    
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/bin/Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/bin/Debug)    
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/bin/Release)  

# cmake module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# C++ version
set (CMAKE_CXX_STANDARD 11)
# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# C++ library postfix
set (CMAKE_DEBUG_POSTFIX "_d")

# OpenCV
if ( WIN32 )
	set(OpenCV_DIR "E:/libs/OpenCV/x64/vc14/lib")
	set(OpenCV_CONFIG_PATH "E:/libs/OpenCV")
	set(WIN_HEADER_PATH "${PROJECT_SOURCE_DIR}/winheadfile")
    include_directories(${WIN_HEADER_PATH})
endif()
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

# Eigen
include_directories("E:/libs/Eigen/include/eigen3")

#cuda
find_package(CUDA)
set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS}
#    -O3 
	-gencode=arch=compute_61,code=\"sm_61,compute_61\"
    )

cuda_add_library(${name} SHARED
    ICCG.cpp
    InputImage.cpp
    BilateralGrid.cpp
)

target_link_libraries(${name}
    ${OpenCV_LIBS}
)

add_executable(bilater BilateralSolver.cpp BilateralSolver.hpp)
target_link_libraries(bilater ${OpenCV_LIBS} )

add_executable(colorize new_color.cpp)
target_link_libraries(colorize ${OpenCV_LIBS} )

add_executable(colorize_test color.cpp)
target_link_libraries(colorize_test ${OpenCV_LIBS})