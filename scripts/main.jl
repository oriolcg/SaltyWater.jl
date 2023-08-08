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

params = NavierStokes_ConvectionDiffusion_params(L=5.0,ne=40,K=1e-6,tf=10.0)
solve_NSCD(params)

end
