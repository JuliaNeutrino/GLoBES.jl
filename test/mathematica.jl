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

    # Blocks: Head[...]{ ... } and Sym{ ... } without a space before {
    b = only(
        parse("MixingMatrix[R12 . R23]{ R12 = Diag[{1, 1, 1}]\n R23 = Rotation[2, 3, theta23] }"),
    )
    @test b.head === :Block && str(b.args[1]) == "MixingMatrix[Dot[R12, R23]]"
    @test length(b.args) == 3 && str(b.args[3]) == "Set[R23, Rotation[2, 3, theta23]]"
    b = only(parse("LaTeX[\"out.tex\"]{ X = Off Print[\"hello\"] Y = On }"))
    @test b.head === :Block && length(b.args) == 4
    @test str(b.args[2]) == "Set[X, Off]" &&
          str(b.args[3]) == "Print[\"hello\"]" &&
          str(b.args[4]) == "Set[Y, On]"
    b = only(parse("Hamiltonian{ H = A }"))
    @test b.head === :Block && b.args[1] === :Hamiltonian
    # With a space, { } is a factor of a product
    @test str(only(parse("ZW = 1/Sqrt[2] {{1, 1}, {I, -I}}"))) ==
          "Set[ZW, Times[Power[Sqrt[2], -1], {{1, 1}, {I, Times[-1, I]}}]]"
    @test only(parse("eta_etau")) === :eta_etau
end
