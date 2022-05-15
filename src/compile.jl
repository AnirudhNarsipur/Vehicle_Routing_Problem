import Pkg
using PackageCompiler

Pkg.add(["JuMP","HiGHS","StatsBase"])
PackageCompiler.create_sysimage([:JuMP,:HiGHS,:StatsBase],
sysimage_path="LSearch.so",
precompile_execution_file="src/precompile.jl")
