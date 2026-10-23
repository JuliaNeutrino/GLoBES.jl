# The help of every exported function shows how to call it, one line per call

@testset "help examples" begin
    meta = Base.Docs.meta(GLoBES)
    for n in names(GLoBES)
        s = String(n)
        (n === :GLoBES || startswith(s, "glb") || startswith(s, "GLB")) && continue
        b = Base.Docs.Binding(GLoBES, n)
        @test haskey(meta, b) || n
        haskey(meta, b) || continue
        doc = join(join(d.text) for d in values(meta[b].docs))
        blocks = [m.captures[1] for m in eachmatch(r"```julia\n(.*?)```"s, doc)]
        @test !isempty(blocks) || n
        for line in split(join(blocks), '\n')
            isempty(strip(line)) && continue
            ex = Meta.parse(line; raise = false)
            @test !(ex isa Expr && ex.head in (:error, :incomplete)) || (n, line)
        end
    end
end
