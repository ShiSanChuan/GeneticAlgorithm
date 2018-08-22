FIND_PACKAGE( OpenCV REQUIRED ) 
FIND_PACKAGE ( Threads REQUIRED )

aux_source_directory(. DIR_SRCS)
add_subdirectory(./cvplotlib)
add_executable(GA ${DIR_SRCS})

TARGET_LINK_LIBRARIES(GA ${OpenCV_LIBS} ${CMAKE_THREAD_LIBS_INIT} cvplotlib)