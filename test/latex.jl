# LaTeX output

@testset "LaTeX" begin
    tex(s) = latex_expression(only(GLoBES.parse_mathematica(s)))
    @test tex("theta12") == "\\theta_{12}"
    @test tex("eta_etau") == "\\eta_{e\\tau}"
    @test tex("Sqrt[Abs[x]] y") == "\\sqrt{\\left|x\\right|} y"
    @test tex("a - b/c") == "a - \\frac{b}{c}"
    @test tex("Dagger[U]") == "{U}^\\dagger"
    @test tex("{{1, 2}, {3, 4}}") == "\\begin{pmatrix} 1 & 2 \\\\ 3 & 4 \\end{pmatrix}"
    @test tex("Exp[I d]") == "e^{i d}"
    @test GLoBES.latex_number(2.5e-5) == "2.5 \\times 10^{-5}"
    @test GLoBES.latex_number(0.3088) == "0.3088"

    m = load_model(joinpath(MODELS, "SNSI"))
    s = latex(m)
    @test occursin("\\section{The Scalar NSI with decoherence model}", s)
    @test occursin("\\frac{1}{2E} M M^\\dagger + H_{\\mathrm{Matter}}", s)
    @test occursin("\\eta_{e\\tau}", s) && occursin("NuFIT 6.1", s)
    @test !occursin("Parameters of the", latex(m; Parameters = false))
    @test occursin("\\toprule", latex(nufit("6.1")))
    @test occursin("\\begin{pmatrix}", latex([1 0; 0 1]))

    r = TexReport("Test"; author = "GLoBES.jl")
    add!(r, m; NuFit = false)
    section!(r, "More")
    table!(r, ["a", "b"], [[1, 2.5e-5]])
    figure!(r, "plot.pdf"; caption = "A plot")
    f = write_tex(tempname() * ".tex", r)
    doc = read(f, String)
    @test startswith(doc, "% Written by GLoBES.jl") && occursin("\\end{document}", doc)
    @test occursin("\\includegraphics[width=0.8\\textwidth]{plot.pdf}", doc)
end
