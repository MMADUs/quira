add_executable(matrix_add
    matrix_add.cpp
)
target_link_libraries(matrix_add
    PRIVATE
        quira::quira
)

# bellstate example
add_executable(bellstate
    bellstate.cpp
)
target_link_libraries(bellstate
    PRIVATE
        quira::quira
)

# teleportation example
add_executable(teleportation
    teleportation.cpp
)
target_link_libraries(teleportation
    PRIVATE
        quira::quira
)

# build examples using `./scripts/build.sh`
#
# run examples:
# ./build/examples/bellstate
# ./build/examples/teleportation
