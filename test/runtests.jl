using Test
using XGrid




########################################################
# Test 


function runtests()

    println("A0: Adjacency transpose from matrix")
    A0=[i+3*(j-1) for i=1:3, j=1:3]
    
    A0t=atranspose(A0)
    println("A0:\n", A0)
    println("A0t:\n", A0t)
    
    @test A0t == VariableTargetAdjacency{Int64}([1, 1, 1, 2, 2, 2, 3, 3, 3], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])


    println("A1: Adjacency transpose from variable target adj")
    A1=VariableTargetAdjacency()
    append!(A1, (1,2,3))
    append!(A1, (4,5,6,7))
    append!(A1, (8,9))
    A1t=atranspose(A1)

    println("A1:\n", A1)
    println("A1t:\n", A1t)


    @test A1t == VariableTargetAdjacency{Int64}([1, 1, 1, 2, 2, 2, 2, 3, 3], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    println("A2: Adjacency transpose from variable target adj")
    A2=VariableTargetAdjacency()
    append!(A2, (1,2,3))
    append!(A2, (2,3,4))
    append!(A2, (3,4,5))
    
    A2t=atranspose(A2)
    println("A2:\n", A2)
    println("A2t:\n", A2t)
    @test A2t == VariableTargetAdjacency{Int64}([1, 1, 2, 1, 2, 3, 2, 3, 3], [1, 2, 4, 7, 9, 10])
end
