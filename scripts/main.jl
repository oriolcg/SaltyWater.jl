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

params = NavierStokes_ConvectionDiffusion_params(
  H=4.0e-2,
  L=1.0,
  μ=1.0e-3,
  𝒟=1.29e-9,
  K=3.6e-3,
  order=2,
  nex=40,ney=20,tf=1.0)
solve_NSCD(params)

end
