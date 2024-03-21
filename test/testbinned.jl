    function testbinning(a)
        dim = size(a, 1)
        n = size(a, 2)
        idx = rand(1:n, n ÷ 4)
        a1 = a[:, idx]

        bpl = BinnedPointList(dim)
        for i = 1:size(a, 2)
            insert!(bpl, a[:, i])
        end

        for i = 1:size(a1, 2)
            ix = insert!(bpl, a1[:, i])
            if ix != idx[i]
                return false
            end
        end
        true
    end
    function testbinning2()
        A = 50
        B = 100
        Z = 20
        d = 5
        builder = SimplexGridBuilder(; Generator = TetGen)
        p1 = point!(builder, 0, 0, 0)
        p2 = point!(builder, B, 0, 0)
        p3 = point!(builder, B, A, 0)
        p4 = point!(builder, B, A + d, 0)
        p5 = point!(builder, B, 2 * A + d, 0)
        p6 = point!(builder, 0, 2 * A + d, 0)
        p7 = point!(builder, 0, A + d, 0)
        p8 = point!(builder, 0, A, 0)
    end

    @test testbinning(rand(1, 10))
    @test testbinning(rand(1, 10000))
    @test testbinning(rand(2, 10))
    @test testbinning(rand(2, 10000))
    @test testbinning(rand(3, 10))
    @test testbinning(rand(3, 10000))

    @test testbinning2() == 8
