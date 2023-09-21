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
  H=7.4e-4,#4.0e-2,
  L=1.5e-2,#1.0,
  μ=8.9e-4,#1.0e-3,
  𝒟=1.5e-9,#1.29e-9,
  K=1/(8.41e10),#3.6e-12,
  ρw=1.0272e3,#1.0e3,
  order=2,
  nex=40,ney=20,tf=2e-5,Δt=1e-5)
solve_NSCD(params)

end
