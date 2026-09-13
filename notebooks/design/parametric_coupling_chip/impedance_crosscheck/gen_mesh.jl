# thin CLI wrapper around mesh_si_cpw_line.jl::generate
#   julia --project=<env-with-Gmsh> gen_mesh.jl OUT.msh [refinement] [length_um] [trace_width_um] [gap_width_um] [ground_width_um]
include(joinpath(@__DIR__, "mesh_si_cpw_line.jl"))

out  = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "stage2_out", "mesh.msh")
refn = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
len  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 4000.0
w    = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 20.0
s    = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 11.5
gw   = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 250.0

mkpath(dirname(out))
generate(; filename=out, order=1, refinement=refn, length_um=len, trace_width_um=w, gap_width_um=s, ground_width_um=gw)
println("wrote ", out)
