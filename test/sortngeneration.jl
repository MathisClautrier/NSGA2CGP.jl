include("../src/NSGA2CGP.jl")
using .NSGA2CGP
using Cambrian
using ArcadeLearningEnvironment
using ArgParse
import Cambrian.mutate
using Test
using Random
import StatsBase.sample

function symbolic_evaluate(i::NSGA2CGPInd; seed::Int=0)
    Random.seed!(seed)
    inputs = rand(i.n_in)
    output = process(i, inputs)
    target = rosenbrock(inputs)
    [-(output[1] - target)^2]
end

function create()
    liste=[]
    n=1
    for i in 1:4
        push!(liste,deepcopy([n,n,n]))
        push!(liste,deepcopy([n+1,n,n]))
        push!(liste,deepcopy([n,n+1,n]))
        push!(liste,deepcopy([n,n,n+1]))
        n+=1
    end
    liste
end

Fitnesses=create()
Ranks=[1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8]

function setfitnesses(e::NSGA2CGPEvolution)
    for i in 1:2*e.config.n_population
        e.population[i].fitness=Fitnesses[i]
    end
end
function selection(pop::Array{NSGA2CGPInd}; t_size=3)
    inds = shuffle!(collect(1:length(pop)))
    sort(pop[inds[1:t_size]])[end]
end

function NSGA2CGPopulate(e::NSGA2CGPEvolution)
    Qt=Array{NSGA2CGPInd}(undef,0)
    for ind in e.population
        push!(Qt,ind)
    end
    for i in 1:e.config.n_population
        p1=selection(e.population)
        child= copy(p1)
        @assert typeof(child)<:NSGA2CGPInd
        #if e.config.p_crossover > 0 && rand() < e.config.p_crossover
        #    parents = vcat(p1, [selection(e.population) for i in 2:e.config.n_parents])
        #    child = crossover(parents...)
        #end
        push!(Qt,child)
    end
    @assert length(Qt)==2*e.config.n_population
    e.population=Qt
end

function fastNonDominatedSort!(e::NSGA2CGPEvolution)

    Fi=Array{NSGA2CGPInd}(undef,0)

    for ind in e.population
        ind.n=0
        ind.r=0
        ind.distance=0.
        ind.S=Array{NSGA2CGPInd}(undef,0)
    end

    for ind1 in e.population
        for ind2 in e.population
            if dominates(e,ind1,ind2)
                push!(ind1.S,ind2)
            elseif dominates(e,ind2,ind1)
                ind1.n+=1
            end
        end
        if ind1.n==0
            ind1.r=1
            push!(Fi,ind1)
        end
    end

    i=1

    while isempty(Fi)==false
        Q=Array{NSGA2CGPInd}(undef,0)
        for ind1 in Fi
            for ind2 in ind1.S
                ind2.n-=1
                if ind2.n==0
                    ind2.r=i+1
                    push!(Q,ind2)
                end
            end
        end
        i=i+1
        Fi=Q
    end
end

function crowdingDistanceAssignement!(e::NSGA2CGPEvolution,I::Array{NSGA2CGPInd})
    for ind in I
        ind.distance=0
    end
    for i in 1:e.config.d_fitness
        sort!(I,by=x->x.fitness[i])
        if I[1].fitness[i]!=I[end].fitness[i]
            I[1].distance=Inf
            I[end].distance=Inf
            quot=I[end].fitness[i]-I[1].fitness[i]
            for j in 2:l-1
                I[j].distance=I[j].distance+
                (I[j+1].fitness[i]-I[j-1].fitness[i])/quot
            end
        end
    end
end

function NSGA2CGPGeneration(e::NSGA2CGPEvolution)
    if(e.gen>1)
        fastNonDominatedSort!(e)
        Pt1=Array{NSGA2CGPInd}(undef,0)
        i=1
        sort!(e.population,by= x -> x.r)
        rank=1
        indIni=1
        indNext=findlast(x -> x.r ==rank,e.population)
        Pt1=[Pt1...,e.population[indIni:indNext]...]
        while indNext < e.config.n_population
            rank+=1
            indIni=indNext
            indNext=findlast(x -> x.r == rank,e.population)
            Pt1=[Pt1...,e.population[indIni+1:indNext]...]
        end
        if isempty(Pt1)
            I=e.population[1:indNext]
            crowdingDistanceAssignement!(e,I)
            sort!(I, by= x->x.distance,rev=true)
            Pt1=I[1:e.config.n_population]
        else
            I=e.population[indIni+1:indNext]
            crowdingDistanceAssignement!(e,I)
            sort!(I, by= x->x.distance,rev=true)
            Pt1=[Pt1...,I[1:e.config.n_population-indIni-1]...]
        end
        @assert length(Pt1)==e.config.n_population
        e.population=Pt1
    end
end

function dominates(e::NSGA2CGPEvolution,ind1::NSGA2CGPInd,ind2::NSGA2CGPInd)
    dom=false
    for i in 1:e.config.d_fitness
        if ind1.fitness[i]<ind2.fitness[i]
            return false
        elseif ind1.fitness[i]>ind2.fitness[i]
            dom=true
        end
    end
    return dom
end

@testset "Test NSGA2 selection and populate" begin
    cfg = get_config("cfg/atari_ram.yaml"; game="test",n_in=128,n_out=18)
    e = NSGA2CGPEvolution(cfg, symbolic_evaluate)
    e.gen=2
    NSGA2CGPopulate(e)
    setfitnesses(e)
    @test length(Ranks)==16
    @test length(e.population)==16
    fastNonDominatedSort!(e)
    sort!(e.population,by= x -> x.r)
    for i in 1:16
        @test Ranks[i]==e.population[i].r
    end
    NSGA2CGPGeneration(e)
    @test length(e.population)==8
    sort!(e.population,by= x -> x.r)
    for i in 1:8
        @test Ranks[i]==e.population[i].r
    end
end
