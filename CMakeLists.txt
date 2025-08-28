cmake_minimum_required(VERSION 3.28)
project(qsfcc)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src/snappy)
add_subdirectory(src/sim_piece)
add_subdirectory(src/sz2)

add_executable(qsfcc
        main.cpp
        src/machete/machete.cpp
        src/machete/huffman.cpp
        src/machete/ovlq.cpp
        src/machete/predict.cpp
        src/machete/hybrid.cpp
        src/qsfcc/qsfcc.h
        src/qsfcc/qsfcc.cpp
        src/qsfcc/sfc.cpp
        src/qsfcc/defs.h
        src/qsfcc/sfc.h
        src/qsfcc/quantizer.h
)

target_include_directories(qsfcc PRIVATE
        ${PROJECT_BINARY_DIR}
        ${PROJECT_BINARY_DIR}/src/snappy
        ${PROJECT_SOURCE_DIR}/src/snappy
        ${CMAKE_SOURCE_DIR}/libs/zstd/include
)

target_link_directories(qsfcc PRIVATE ${CMAKE_SOURCE_DIR}/libs/zstd/dll)

target_link_libraries(qsfcc PRIVATE snappy zstd sim_piece sz)