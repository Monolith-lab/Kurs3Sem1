cmake_minimum_required(VERSION 4.0)
project(Raschetnoe_zadanie_kotli C)

set(CMAKE_C_STANDARD 11)

add_executable(Raschetnoe_zadanie_kotli main.c)

target_link_libraries(Raschetnoe_zadanie_kotli PRIVATE
    m)