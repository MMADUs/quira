cmake_minimum_required(VERSION 3.20)

project(
    quira
    VERSION 0.1.0
    DESCRIPTION "General purpose quantum computing package"
    LANGUAGES CXX
)

# using C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# build examples & tests, enabled by default
# cmake -S . -B build -D BUILD_EXAMPLES=OFF
option(QUIRA_BUILD_EXAMPLES "Build example programs" ON)
option(QUIRA_BUILD_TESTS "Build test programs" ON)

# get eigen dependency
find_package(Eigen3 CONFIG QUIET)

if(NOT Eigen3_FOUND)
    include(FetchContent)

    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )

    FetchContent_GetProperties(eigen)

    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()

    add_library(Eigen3::Eigen INTERFACE IMPORTED)

    target_include_directories(Eigen3::Eigen
        INTERFACE
            ${eigen_SOURCE_DIR}
    )
endif()

# collect all .cpp files inside src/*
file(GLOB_RECURSE QUIRA_SOURCES CONFIGURE_DEPENDS
    src/*.cpp
)
# main library target
add_library(quira
    ${QUIRA_SOURCES}
)

# namespace alias
add_library(quira::quira ALIAS quira)

# target header files
target_include_directories(quira
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# link eigen to quira
target_link_libraries(quira
    PUBLIC
        Eigen3::Eigen
)

# require C++20 for this target
target_compile_features(quira
    PUBLIC
        cxx_std_20
)

if(QUIRA_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(QUIRA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
