# The help of every exported function shows how to call it, one line per call

using Markdown

# The julia code blocks of a docstring. A docstring copied with
# `@doc (@doc f) g` keeps only its Markdown on Julia 1.10, not its text.
function julia_blocks(d::Base.Docs.DocStr)
    md = d.object isa Markdown.MD ? d.object : Markdown.parse(join(string.(d.text)))
    blocks = String[]
    walk(x) = nothing
    walk(x::Markdown.MD) = foreach(walk, x.content)
    walk(x::Markdown.Code) = x.language == "julia" && push!(blocks, x.code)
    walk(md)
    return blocks
end

@testset "help examples" begin
    meta = Base.Docs.meta(GLoBES)
    exported = filter(names(GLoBES)) do n
        s = String(n)
        !(n === :GLoBES || startswith(s, "glb") || startswith(s, "GLB"))
    end
    bindings = [Base.Docs.Binding(GLoBES, n) for n in exported]
    @test [n for (n, b) in zip(exported, bindings) if !haskey(meta, b)] == Symbol[]   # no help
    documented = [(n, b) for (n, b) in zip(exported, bindings) if haskey(meta, b)]
    examples = Dict(n => reduce(vcat, julia_blocks.(values(meta[b].docs))) for (n, b) in documented)
    @test [n for (n, _) in documented if isempty(examples[n])] == Symbol[]   # help without examples
    # Each line is one complete call
    broken = [
        (n, line) for (n, blocks) in examples for
        line in split(join(blocks, '\n'), '\n') if !isempty(strip(line)) && (
            ex = Meta.parse(line; raise = false);
            ex isa Expr && ex.head in (:error, :incomplete)
        )
    ]
    @test broken == Tuple{Symbol, SubString{String}}[]
end
