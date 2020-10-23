using Test
include("../src/NSGA2CGP.jl")
using .NSGA2CGP
import YAML

@testset "CGPInd construction" begin
    cfg = get_config("test/test.yaml")
    ind = NSGA2CGPInd(cfg)

    @test length(ind.nodes) == 3 * 10 + 4
    for node in ind.nodes
        if node.active
            @test node.x >= 1
            @test node.x <= length(ind.nodes)
            @test node.y >= 1
            @test node.y <= length(ind.nodes)
        end
        # test that stringifying works
        @test typeof(string(node)) == String
    end
end

"""
using random values, sort individuals that are different
"""
function select_random(pop::Array{NSGA2CGPInd}, elite::Int; n_in=113, n_sample=100)
    actions = zeros(Int, length(pop))
    dists = zeros(n_sample, length(pop))
    inputs = rand(n_in, n_sample)

    for i in 1:n_sample
        for j in eachindex(pop)
            actions[j] = argmax(process(pop[j], inputs[:, i]))
        end
        for j in eachindex(pop)
            dists[i, j] = sum(actions[j] .!= actions)
        end
    end
    d = sum(dists, dims=1)[:]
    ds = sortperm(d)[1:elite]
    pop[ds]
end

@testset "Processing" begin
    cfg = get_config("test/test.yaml"; functions=["f_abs", "f_add", "f_mult"])
    ind = NSGA2CGPInd(cfg)

    # test that f(0, 0, 0, 0) = 0
    inputs = zeros(4)
    set_inputs(ind, inputs)
    for i in 1:4
        @test ind.buffer[i] == 0.0
    end
    output = process(ind)
    @test output[1] == 0.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] == 0.0
        end
    end

    # test that f(1, 1, 1, 1) = 1
    for i in eachindex(ind.nodes)
        ind.buffer[i] = 1.0 # requires that buffer is 1
    end
    output = process(ind, ones(4))
    @test output[1] == 1.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] == 1.0
        end
    end

    cfg = get_config("test/test.yaml"; functions=["f_abs", "f_add", "f_mult"], recur=1.0)
    ind = NSGA2CGPInd(cfg)
    output = process(ind, rand(4))
    @test output[1] <= 1.0 && output[1] >= -1.0
    for i in eachindex(ind.nodes)
        if ind.nodes[i].active
            @test ind.buffer[i] <= 1.0 && ind.buffer[i] >= -1.0
        end
    end

    pop = [NSGA2CGPInd(cfg) for i in 1:10]
    sp = select_random(pop, 2; n_in=cfg.n_in, n_sample=5)
    @test length(sp) == 2
    @test sp[1].buffer[1] != 0.0
end
