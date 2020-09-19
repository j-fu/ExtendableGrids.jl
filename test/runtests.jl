using Test


using Test

modname(fname)=splitext(basename(fname))[1]

#
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method test() returning true
# or false depending on success.
#
function run_tests_from_directory(testdir,prefix)
    println("Directory $(testdir):")
    begin
        examples=modname.(readdir(testdir))
        for example in examples
            if length(example)>=length(prefix) &&example[1:length(prefix)]==prefix
                @info "$(example):"
                path=joinpath(testdir,"$(example).jl")
                @eval begin
                    include($path)
                    # Compile + run test
                    @info "compile:"
                    @time @test eval(Meta.parse("$($example).test()"))
                    # Second run: pure execution time.
                    @info "run:"
                    @time eval(Meta.parse("$($example).test()"))
                end
            end
        end
    end
end


function run_all_tests()
    @time begin
        run_tests_from_directory(@__DIR__,"test_")
        run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example")
        @info "Tests finished"
    end
end

run_all_tests()

# include("adjacency.jl")
# include("extendablegrid.jl")


# test_adj_correctness()
# xadj=test_lazy(test_create_square())

# @test xadj[:,5]==Int32[46, 80, 88, 118, 159, 162, 165, 167]

