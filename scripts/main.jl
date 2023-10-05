module Main

using DrWatson
@quickactivate "SaltyWater.jl"

# Here you may include files from the source directory
include(srcdir("NavierStokes_ConvectionDiffusion.jl"))
include(srcdir("NavierStokes_ConvectionDiffusion_static.jl"))
include(srcdir("NavierStokes_ConvectionDiffusion_static_withWT.jl"))
using .NavierStokes_ConvectionDiffusion
using .NavierStokes_ConvectionDiffusion_Static
using .NavierStokes_ConvectionDiffusion_static_withWT

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)

# params = NavierStokes_ConvectionDiffusion_params()
# solve_NSCD(params)

# params = NavierStokes_ConvectionDiffusion_params(
#   H=7.4e-4,#4.0e-2,
#   L=1.5e-2,#1.0,
#   μ=8.9e-4,#1.0e-3,
#   𝒟=1.5e-9,#1.29e-9,
#   K=1/(8.41e10),#3.6e-12,
#   ρw=1.0272e3,#1.0e3,
#   order=2,
#   nex=40,ney=40,
#   ϕ∞=35000.0)
# solve_NSCD(params)

#params = NavierStokes_ConvectionDiffusion_static_params(ney=6)
#solve_NSCD_static(params)

params = NavierStokes_ConvectionDiffusion_static_withWT_params(ney=6)
solve_NSCD_static_withWT(params)

# params = NavierStokes_ConvectionDiffusion_static_params(
#   H=7.4e-4,
#   L=1.5e-2,
#   μ=8.9e-4,
#   ρw=1.027e3,
#   𝒟=1.5e-9,
#   U∞=0.129,
#   ϕ∞=600,
#   order=2,
#   nex=200,ney=40,
# )
# solve_NSCD_static(params)

params = NavierStokes_ConvectionDiffusion_static_withWT_params(
  H=7.4e-4,
  L=1.5e-2,
  μ=8.9e-4,
  ρw=1.027e3,
  𝒟=1.5e-9,
  U∞₀=0.129,
  ϕ∞=600,
  order=2,
  nex=100,ney=20,
  Δt = 1.0e-2,
  tf=1.0e-1
)
solve_NSCD_static_withWT(params)

end
