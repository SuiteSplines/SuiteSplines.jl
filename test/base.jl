using Test, SafeTestsets

@safetestset "Even or odd numbers" begin

    using IgaBase

    @test select_even_or_odd(4) == Even()
    @test select_even_or_odd(5) == Odd()

    @test_throws NotImplementedError throw(NotImplementedError("he"))
end
