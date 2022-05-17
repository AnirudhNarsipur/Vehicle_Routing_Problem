import Pkg; Pkg.add("PackageCompiler")
using PackageCompiler

Pkg.add(["JuMP","HiGHS","StatsBase"])
PackageCompiler.create_sysimage([:JuMP,:HiGHS,:StatsBase],
sysimage_path="src/LSearch.so",
precompile_execution_file="src/precompile.jl")
