include("../src/NSGA2CGP.jl")
using .NSGA2CGP
using Cambrian
using ArcadeLearningEnvironment
using ArgParse
import Cambrian.mutate
import Random


s = ArgParseSettings()
@add_arg_table s begin
    "--cfg"
    help = "configuration script"
    default = "cfg/atari_ram.yaml"
    "--game"
    help = "game rom name"
    arg_type = Array{String}
    default = ["centipede","frostbite","ms_pacman"]
    "--seed"
    help = "random seed for evolution"
    arg_type = Int
    default = 0
    "--ind"
    help = "individual for evaluation"
    arg_type = String
    default = ""
end
args = parse_args(ARGS, s)

cfg = get_config(args["cfg"]; game=args["game"], n_in=3, n_out=3)
Random.seed!(args["seed"])

mutate(i::NSGA2CGPInd) = goldman_mutate(cfg, i)
fit(i::NSGA2CGPInd) = 10
e = NSGA2CGPEvolution(cfg, fit)
