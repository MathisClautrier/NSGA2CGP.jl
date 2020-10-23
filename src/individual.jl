export Node, NSGA2CGPInd
import Base.copy, Base.String, Base.show, Base.summary
import Cambrian.print

#Been modified
"default function for nodes, will cause error if used as a function node"
function f_null(args...)::Nothing
    nothing
end

struct Node
    x::Int16
    y::Int16
    f::Function
    active::Bool
end

mutable struct NSGA2CGPInd <: Cambrian.Individual
    n_in::Int16
    n_out::Int16
    chromosome::Array{Float64}
    genes::Array{Int16}
    outputs::Array{Int16}
    nodes::Array{Node}
    buffer::Array{Float64}
    fitness::Array{Float64}
    distance::Float64
    r::Int64
    n::Int64
    S::Array{NSGA2CGPInd}
end

function recur_active!(active::BitArray, ind::Int16, xs::Array{Int16},
                       ys::Array{Int16}, fs::Array{Int16},
                       two_arity::BitArray)::Nothing
    if ind > 0 && ~active[ind]
        active[ind] = true
        recur_active!(active, xs[ind], xs, ys, fs, two_arity)
        if two_arity[fs[ind]]
            recur_active!(active, ys[ind], xs, ys, fs, two_arity)
        end
    end
end

function find_active(cfg::NamedTuple, genes::Array{Int16},
                     outputs::Array{Int16})::BitArray
    R = cfg.rows
    C = cfg.columns
    active = falses(R, C)
    xs = genes[:, :, 1] .- Int16(cfg.n_in)
    ys = genes[:, :, 2] .- Int16(cfg.n_in)
    fs = genes[:, :, 3]
    for i in eachindex(outputs)
        recur_active!(active, outputs[i] - Int16(cfg.n_in), xs, ys, fs,
                      cfg.two_arity)
    end
    active
end

function NSGA2CGPInd(cfg::NamedTuple, chromosome::Array{Float64}, genes::Array{Int16},
                outputs::Array{Int16})::NSGA2CGPInd
    R = cfg.rows
    C = cfg.columns
    nodes = Array{Node}(undef, R * C + cfg.n_in)
    for i in 1:cfg.n_in
        nodes[i] = Node(0, 0, f_null, false)
    end
    i = cfg.n_in
    active = find_active(cfg, genes, outputs)
    for y in 1:C
        for x in 1:R
            i += 1
            nodes[i] = Node(genes[x, y, 1], genes[x, y, 2],
                            cfg.functions[genes[x, y, 3]],
                            active[x, y])
        end
    end
    buffer = zeros(R * C + cfg.n_in)
    fitness = -Inf .* ones(cfg.d_fitness)
    NSGA2CGPInd(cfg.n_in, cfg.n_out, chromosome, genes, outputs, nodes, buffer, fitness,0.,0,0,[])
end

function NSGA2CGPInd(cfg::NamedTuple, chromosome::Array{Float64})::NSGA2CGPInd
    R = cfg.rows
    C = cfg.columns
    genes = reshape(chromosome[1:(R*C*3)], R, C, 3)
    # TODO: recurrency is ugly and slow
    maxs = collect(1:R:R*C)
    maxs = round.((R*C .- maxs) .* cfg.recur .+ maxs)
    maxs = min.(R*C + cfg.n_in, maxs .+ cfg.n_in)
    maxs = repeat(maxs, 1, R)'
    genes[:, :, 1] .*= maxs
    genes[:, :, 2] .*= maxs
    genes[:, :, 3] .*= length(cfg.functions)
    genes = Int16.(ceil.(genes))
    outputs = Int16.(ceil.(chromosome[(R*C*3+1):end] .* (R * C + cfg.n_in)))
    NSGA2CGPInd(cfg, chromosome, genes, outputs)
end

function NSGA2CGPInd(cfg::NamedTuple)::NSGA2CGPInd
    chromosome = rand(cfg.rows * cfg.columns * 3 + cfg.n_out)
    NSGA2CGPInd(cfg, chromosome)
end

function NSGA2CGPInd(cfg::NamedTuple, ind::String)::NSGA2CGPInd
    dict = JSON.parse(ind)
    NSGA2CGPInd(cfg, Array{Float64}(dict["chromosome"]))
end

function copy(n::Node)
    Node(n.x, n.y, n.f, n.active)
end

function copy(ind::NSGA2CGPInd)
    buffer = Array{Float64}(undef, length(ind.buffer))
    nodes = Array{Node}(undef, length(ind.nodes))
    for i in eachindex(ind.nodes)
        nodes[i] = copy(ind.nodes[i])
    end
    NSGA2CGPInd(ind.n_in, ind.n_out, copy(ind.chromosome), copy(ind.genes),
           copy(ind.outputs), nodes, buffer, copy(ind.fitness),copy(ind.distance),
           copy(ind.r),copy(ind.n),copy(ind.S))
end

function String(n::Node)
    JSON.json(Dict(:x=>n.x, :y=>n.y, :f=>string(n.f), :active=>n.active))
end

function show(io::IO, n::Node)
    print(io, String(n))
end

function String(ind::NSGA2CGPInd)
    JSON.json(Dict("chromosome"=>ind.chromosome, "fitness"=>ind.fitness))
end

function show(io::IO, ind::NSGA2CGPInd)
    print(io, String(ind))
end

function get_active_nodes(ind::NSGA2CGPInd)
    ind.nodes[[n.active for n in ind.nodes]]
end

function rank!(ind::NSGA2CGPInd,rank::Int64)
    ind.r=rank
end

function n_increment!(ind::NSGA2CGPInd,m::Int64)
    ind.n+=m
end

function distance!(ind::NSGA2CGPInd,dist::Float64)
    ind.distance=dist
end

function addDominated!(ind::NSGA2CGPInd,dominated::NSGA2CGPInd)
    push!(ind.S,dominated)
end

function reset_params!(ind::NSGA2CGPInd)
    ind.n=0
    ind.r=0
    ind.distance=0
    ind.S=Array{NSGA2CGPInd}(undef,0)
end

function get_S(ind::NSGA2CGPInd)
    ind.S
end

function get_n(ind::NSGA2CGPInd)
    ind.n
end

function get_distance(ind::NSGA2CGPInd)
    ind.distance
end

function get_rank(ind::NSGA2CGPInd)
    ind.r
end

function get_fitness_i(ind::NSGA2CGPInd,i::Int64)
    ind.fitness[i]
end


function summary(io::IO, ind::NSGA2CGPInd)
    print(io, string("NSGA2CGPInd(", get_active_nodes(ind), ", ",
                     findall([n.active for n in ind.nodes]), ", ",
                     ind.outputs, " ,",
                     ind.fitness, ")"))
end

function interpret(i::NSGA2CGPInd)
    x::AbstractArray->process(i, x)
end
