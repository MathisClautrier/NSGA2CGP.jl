export NSGA2CGPEvolution

import Cambrian.populate, Cambrian.evaluate, Cambrian.selection, Cambrian.generation
import StatsBase.sample

#been modified
mutable struct NSGA2CGPEvolution{T} <: Cambrian.AbstractEvolution
    config::NamedTuple
    logger::CambrianLogger
    population::Array{T}
    fitness::Function
    gen::Int
end

function tournamentSelection(pop::Array{NSGA2CGPInd}; t_size=3)
    inds = shuffle!(collect(1:length(pop)))
    sort(pop[inds[1:t_size]])[end]
end

function NSGA2CGPPopulate(e::NSGA2CGPEvolution)
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
        if e.config.p_mutation > 0 && rand() < e.config.p_mutation
            child = mutate(child)
        end
        push!(Qt,child)
    end
    @assert length(Qt)==2*e.config.n_population
    e.population=Qt
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

populate(e::NSGA2CGPEvolution) = NSGA2CGPPopulate(e)
evaluate(e::NSGA2CGPEvolution) = Cambrian.fitness_evaluate(e, e.fitness)
generation(e::NSGA2CGPEvolution) = NSGA2CGPGeneration(e)

function NSGA2CGPEvolution(cfg::NamedTuple, fitness::Function;
                      logfile=string("logs/", cfg.id, ".csv"))
    logger = CambrianLogger(logfile)
    population = Cambrian.initialize(NSGA2CGPInd, cfg)
    NSGA2CGPEvolution(cfg, logger, population, fitness, 0)
end
