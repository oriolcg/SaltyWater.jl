module Main

using DrWatson
@quickactivate "SaltyWater.jl"

# Here you may include files from the source directory
include(srcdir("NavierStokes_ConvectionDiffusion.jl"))
using .NavierStokes_ConvectionDiffusion

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)

params = NavierStokes_ConvectionDiffusion_params()
solve_NSCD(params)

params = NavierStokes_ConvectionDiffusion_params(nex=40,ney=80,tf=1.0,μ=1.0e-3)
solve_NSCD(params)

end
