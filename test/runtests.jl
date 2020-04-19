using Test
using XGrid

include("adjacency.jl")
include("extendablegrid.jl")


test_adj_correctness()
xadj=test_lazy(test_create_square())

@test xadj[:,5]==Int32[46, 80, 88, 118, 159, 162, 165, 167]

