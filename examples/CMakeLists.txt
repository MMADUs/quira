add_executable(matrix_add
    matrix_add.cpp
)

target_link_libraries(matrix_add
    PRIVATE
        quira::quira
)
