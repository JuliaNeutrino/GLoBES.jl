# The parser of Mathematica-like model files

@testset "Mathematica parser" begin
    parse = GLoBES.parse_mathematica
    MExpr = GLoBES.MExpr
    str(x) = sprint(GLoBES.show_mexpr, x)

    @test str(only(parse("a + b c^2 - d/e"))) ==
          "Plus[a, Times[b, Power[c, 2]], Times[-1, d, Power[e, -1]]]"
    @test str(only(parse("f[x, {1, 2}] := x . y"))) == "SetDelayed[f[x, {1, 2}], Dot[x, y]]"
    @test only(parse("1.5*^-3")) == 1.5e-3
    @test str(only(parse("Model`Name = \"X\" (* a (* nested *) comment *)"))) ==
          "Set[Model`Name, \"X\"]"
    @test length(parse("a = 1\nb = 2;\nc = {1,\n 2}")) == 3
    @test str(only(parse("x \\[RightArrow] y"))) == "Rule[x, y]"

    # SARAH-style entries
    @test str(only(parse("AdditionalTerms[[1]] = {Matter, A Diag[{1, 0, 0}]};"))) ==
          "Set[Part[AdditionalTerms, 1], {Matter, Times[A, Diag[{1, 0, 0}]]}]"
    # { } is a factor of a product
    @test str(only(parse("ZW = 1/Sqrt[2] {{1, 1}, {I, -I}}"))) ==
          "Set[ZW, Times[Power[Sqrt[2], -1], {{1, 1}, {I, Times[-1, I]}}]]"
    @test only(parse("eta_etau")) === :eta_etau
end
