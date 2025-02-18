module Main

using DrWatson
@quickactivate "SaltyWater.jl"

using CSV
using DataFrames

# Here you may include files from the source directory
# include(srcdir("NavierStokes_ConvectionDiffusion.jl"))
# include(srcdir("NavierStokes_ConvectionDiffusion_static.jl"))
include(srcdir("NavierStokes_ConvectionDiffusion_static_pout.jl"))
# include(srcdir("NavierStokes_ConvectionDiffusion_static_withWT.jl"))
#using .NavierStokes_ConvectionDiffusion
#using .NavierStokes_ConvectionDiffusion_Static
using .NavierStokes_ConvectionDiffusion_Static_pout
#using .NavierStokes_ConvectionDiffusion_static_withWT

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

params = NavierStokes_ConvectionDiffusion_static_params(nex=12,ney=3)
_,_,_ = solve_NSCD_static(params)

#params = NavierStokes_ConvectionDiffusion_static_withWT_params(ney=6)
#solve_NSCD_static_withWT(params)

p_out_vec = [8000]
u_vec = Float64[]
ϕ_vec = Float64[]
p_vec = Float64[]

for p_out in p_out_vec
  params = NavierStokes_ConvectionDiffusion_static_params(
    H=7.4e-4,
    L=(1.5e-1)/2,
    μ=8.9e-4,
    ρw=1.027e3,
    𝒟=1.5e-9,
    U∞=  0.09, #0.129, #,0.258, #
    ϕ∞=600,
    order=2,
    nex=600,ney=15,
    pₒ= p_out,
  )
  u, ϕ, p = solve_NSCD_static(params)
  push!(u_vec,u)
  push!(ϕ_vec,ϕ)
  push!(p_vec,p)
end

 file = open("./data.csv","w")
 df = DataFrame(u=u_vec, ϕ = ϕ_vec, p=p_vec, p_out = p_out_vec)
 @show df
 CSV.write(file,df)
 
 close(file)



#params = NavierStokes_ConvectionDiffusion_static_withWT_params(
  #H=7.4e-4,
  #L=1.5e-2,
  #μ=8.9e-4,
  #ρw=1.027e3,
  #𝒟=1.5e-9,
  #U∞₀=0.129,
  #ϕ∞=600,
  #order=2,
  #nex=100,ney=20,
  #Δt = 1.0e-2,
  #tf=1.0e-1
#)
#solve_NSCD_static_withWT(params)

end
