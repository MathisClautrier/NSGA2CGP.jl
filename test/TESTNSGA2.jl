using Test
include("../src/NSGA2CGP.jl")
using .NSGA2CGP
import YAML

@testset "NSGA2CGPInd construction" begin
    cfg = get_config("cfg/atari_ram.yaml"; game="test",n_in=128,n_out=18)
    ind = NSGA2CGPInd(cfg)

    @test length(ind.nodes) == 1 * 100 + 128
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

@testset "New features for NSGA2CGP" begin
    cfg = get_config("cfg/atari_ram.yaml"; game="test",n_in=128,n_out=18)
    ind1 = NSGA2CGPInd(cfg)
    ind2 = copy(ind1)
    @test ind1.n_in==128
    @test ind1.n_out==18
    @test ind1.chromosome==ind2.chromosome
    @test ind1.genes==ind2.genes
    @test ind1.outputs==ind2.outputs
    @test ind1.nodes==ind2.nodes
    @test ind1.fitness==ind2.fitness
    @test ind1.distance==0
    @test ind1.r==0
    @test ind1.S==[]
    @test ind1.n==0
    @test typeof(ind1.S) <: Array{NSGA2CGPInd}
    ind1.n+=1
    @test ind1.n==1
    push!(ind1.S,ind2)
    @test !isempty(ind1.S)
    @test ind1.S[1]==ind2
    ind1.r,ind2.r,ind3.r,ind4.r=1,2,2,3
    pop=[]
end

@testset "New methods for NSGA2CGP" begin
    cfg = get_config("cfg/atari_ram.yaml"; game="test",n_in=128,n_out=18)
    pop=[NSGA2CGPInd(cfg) for i in 1:10]
    ranks=[1,2,2,3,4,5,5,6,7,8]
    for i in range(10):
        pop[i].r=ranks[i]
    end
    @test findlast(x,)
end
