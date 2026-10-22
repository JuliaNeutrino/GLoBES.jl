# Control files (Run.m)

@testset "run_model" begin
    out = mktempdir()
    r = run_model(joinpath(MODELS, "SNSI"); output = out, plots = false)
    @test length(r.calculations) == 4
    @test r.calculations[4][2] ≈ [0.005470290919721274, 0.023351950381050646, 0.010313694447619093] rtol =
        1e-10
    @test r.tex == joinpath(out, "SNSI.tex") && isfile(r.tex)
    doc = read(r.tex, String)
    @test occursin("Effective mixing matrix", doc)
    @test occursin("diagonalised numerically", doc)               # Print[...]
    @test !occursin("NuFIT 6.1, normal ordering", doc)           # NuFit = Off

    # A fit: Δχ² of η_eτ at NOvA, the data simulated without NSI
    dir = mktempdir()
    cp(joinpath(MODELS, "SNSI", "SNSI.m"), joinpath(dir, "SNSI.m"))
    cp(joinpath(MODELS, "SNSI", "parameters.m"), joinpath(dir, "parameters.m"))
    nova = joinpath(DATADIR, "NOvA.glb")
    write(
        joinpath(dir, "Run.m"),
        """
        Fit[Eta]{ Experiment = "$nova", Scan = {eta_etau -> Range[-0.1, 0.1, 0.05]}, Output = "eta.dat" }
        LaTeX["fit.tex"]{ Model = Off, Parameters = Off }
        """,
    )
    r = run_model(dir; plots = false)
    name, res = only(r.fits)
    @test name == "Eta" && res.x ≈ -0.1:0.05:0.1
    @test res.chi2[3] == 0 && all(res.chi2[[1, 2, 4, 5]] .> 0)
    @test length(readlines(joinpath(dir, "eta.dat"))) == 6
    @test occursin("\\paragraph{Eta}", read(joinpath(dir, "fit.tex"), String))
end
