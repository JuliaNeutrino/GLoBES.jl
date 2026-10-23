# Particle models from model files: reading, symmetry checks and the
# tree-level spectrum

const SCOTOGENIC = """
Model`Name = "Scotogenic";
Global[[1]] = {Z[2], Z2};
Gauge[[1]] = {B,  U[1], hypercharge, g1, False, 1};
Gauge[[2]] = {WB, SU[2], left,       g2, True,  1};
Gauge[[3]] = {G,  SU[3], color,      g3, False, 1};
FermionFields[[1]] = {q, 3, {uL, dL}, 1/6, 2, 3, 1};
FermionFields[[2]] = {l, 3, {vL, eL}, -1/2, 2, 1, 1};
FermionFields[[3]] = {d, 3, conj[dR], 1/3, 1, -3, 1};
FermionFields[[4]] = {u, 3, conj[uR], -2/3, 1, -3, 1};
FermionFields[[5]] = {e, 3, conj[eR], 1, 1, 1, 1};
FermionFields[[6]] = {n, 3, conj[nR], 0, 1, 1, -1};
ScalarFields[[1]] = {H, 1, {Hp, H0}, 1/2, 2, 1, 1};
ScalarFields[[2]] = {Et, 1, {etp, et0}, 1/2, 2, 1, -1};
NameOfStates = {GaugeES, EWSB};
DEFINITION[GaugeES][LagrangianInput] = {
    {LagFer, {AddHC->True}}, {LagNV, {AddHC->True}}, {LagH, {AddHC->False}},
    {LagEt, {AddHC->False}}, {LagHEt, {AddHC->False}}, {LagHEtHC, {AddHC->True}}
};
LagFer = Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q + Yn Et.n.l;
LagNV = Mn/2 n.n;
LagH = -(-mH2 conj[H].H + 1/2 lambda1 conj[H].H.conj[H].H);
LagEt = -(+mEt2 conj[Et].Et + 1/2 lambda2 conj[Et].Et.conj[Et].Et);
LagHEt = -(+lambda3 conj[H].H.conj[Et].Et + lambda4 conj[H].Et.conj[Et].H);
LagHEtHC = -(+1/2 lambda5 conj[H].Et.conj[H].Et);
DEFINITION[EWSB][VEVs] = {
    {H0,  {v, 1/Sqrt[2]}, {Ah, I/Sqrt[2]}, {hh, 1/Sqrt[2]}},
    {et0, {0, 0}, {etI, I/Sqrt[2]}, {etR, 1/Sqrt[2]}}
};
"""

function scotogenic_model(extra = "")
    m = BSMModel("Scotogenic")
    GLoBES.read_particle_string!(m, SCOTOGENIC * extra)
    return m
end

const SCOTO_VALUES = Dict{Symbol, Any}(
    :v => 246.22,
    :lambda1 => 0.26,
    :lambda2 => 0.1,
    :lambda3 => 0.2,
    :lambda4 => -0.1,
    :lambda5 => 1e-9,
    :mEt2 => 1e6,
    :Mn => Diagonal([1e3, 2e3, 3e3]),
)

@testset "particle models" begin
    m = scotogenic_model()
    @test m.name == "Scotogenic"
    @test length(m.fields) == 8 && length(m.gauge) == 3
    c = check_model(m)
    @test GLoBES.passed(c)

    # A Yukawa coupling that breaks Z2 and one that breaks hypercharge
    bad = scotogenic_model(
        "LagBad = Yb H.n.l;\nLagBad2 = Yc H.e.l;\nDEFINITION[GaugeES][LagrangianInput] = {{LagBad, {AddHC->True}}, {LagBad2, {AddHC->True}}};",
    )
    c = check_model(bad)
    @test length(c.problems) == 2
    @test any(p -> occursin("Z2", p), c.problems) &&
          any(p -> occursin("hypercharge", p), c.problems)

    # The same model built in Julia
    j = BSMModel("Scoto")
    add_gauge!(j, :B, :U1, :hypercharge, :g1)
    add_gauge!(j, :WB, (:SU, 2), :left, :g2; expand = true)
    add_gauge!(j, :G, (:SU, 3), :color, :g3)
    add_global!(j, :Z2; Z = 2)
    add_fermion!(j, :l, 3, [:vL, :eL], -1 // 2, 2, 1, 1)
    add_fermion!(j, :n, 3, :nR, 0, 1, 1, -1; conjugated = true)
    add_scalar!(j, :H, 1, [:Hp, :H0], 1 // 2, 2, 1, 1)
    add_scalar!(j, :Et, 1, [:etp, :et0], 1 // 2, 2, 1, -1)
    add_lagrangian!(j, :LagNu, "Yn Et.n.l + Mn/2 n.n"; hc = true)
    @test GLoBES.passed(check_model(j))

    # Tree-level spectrum: m_h² = λ1 v², m²(η_R,I) = m²_η + (λ3 + λ4 ± λ5) v²/2
    s = spectrum(m, SCOTO_VALUES; tadpoles = [:mH2])
    v = 246.22
    @test GLoBES.mass(s, :hh) ≈ sqrt(0.26) * v rtol = 1e-6
    @test GLoBES.mass(s, :etR) ≈ sqrt(1e6 + (0.2 - 0.1 + 1e-9) * v^2 / 2) rtol = 1e-6
    @test GLoBES.mass(s, :etI) ≈ sqrt(1e6 + (0.2 - 0.1 - 1e-9) * v^2 / 2) rtol = 1e-6
    @test GLoBES.mass(s, :etp) ≈ sqrt(1e6 + 0.2 * v^2 / 2) rtol = 1e-6
    @test GLoBES.mass(s, :Ah) == 0
    @test s.tadpoles[:mH2] ≈ 0.26 * v^2 / 2 rtol = 1e-6

    # SLHA files
    m = scotogenic_model("ParameterDefinitions = {{Yn, {LesHouches -> YN}}};")
    f = tempname() * ".slha"
    Y = [0.1 0.2 0.0; 0.0 0.3 0.1im; 0.05 0.0 0.4]
    write_leshouches(f, m, Dict(:Yn => Y))
    @test haskey(read_leshouches(f), "YN")
    @test leshouches_values(m, f)[:Yn] ≈ Y
end
