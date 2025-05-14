# pre-setup Cargo workspace

## in-case i forgot how to setup workspaces

add the bin executable
> cargo new <bin-name> --bin --vcs none

add the cargo lib
> cargo new <lib-name> --lib --vcs none

running the workspace from bin
> cargo run -p <bin-name>
