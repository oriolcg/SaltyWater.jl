module NavierStokes_ConvectionDiffusion
"""
    NavierStokes_ConvectionDiffusion

    This module contains the functions for the Navier-Stokes equations with convection-diffusion.

    # Examples
    ```julia-repl
    julia> using NavierStokes_ConvectionDiffusion
    julia> solve_NSCD()
    ```
"""

using Gridap
using Parameters
using DrWatson

export NavierStokes_ConvectionDiffusion_params, solve_NSCD

"""
    solve_NSCD(params)

    This function solves the Navier-Stokes equations with convection-diffusion on
    a 2D rectangular domain using the finite element method. We assume a channel
    flow with a parabolic inflow profile and an osmosis membrane in the middle of
    the channel. The membrane is modeled as a Robin boundary condition.

    # Arguments
    - `params::Parameters`: A Parameters object containing the parameters of the problem.

    # Output
    - vtk file with the solution, stored in the `data/sims` folder.
"""
function solve_NSCD(params)

  # Define the domain
  @unpack H,L,ne = params
  𝒯 = CartesianDiscreteModel((0,L,0,2H), (ne,ne))
  Ω = Interior(𝒯)

  # Divide channel
  function is_in_feed_channel(coords)
    coords_avg = sum(coords)/length(coords) # center of the cell
    [ xᵢ[2] < H for xᵢ in coords_avg] # check if the center of the cell is above or below the membrane
  end
  coords = get_cell_coordinates(Ω)
  feed_channel_mask = is_in_feed_channel(coords)
  feed_channel_indeces = findall(feed_channel_mask)
  permeate_channel_indeces = findall(.!feed_channel_mask)
  Ωf = Interior(Ω, feed_channel_indeces)
  Ωp = Interior(Ω, permeate_channel_indeces)

  # Define boundary tags
  labels_Ω = get_face_labeling(𝒯)
  add_tag_from_tags!(labels_Ω,"top",[3,4,6])       # assign the label "top" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)

  # Define boundaries
  Γfp = Interface(Ωf,Ωp)
  Γin = Boundary(Ω, tags="inlet")
  Γout = Boundary(Ω, tags="outlet")
  Γtop = Boundary(Ω, tags="top")
  Γb = Boundary(Ω, tags="bottom")
  nfp = get_normal_vector(Γfp)

  # Boundary condition
  @unpack U∞ = params
  uin((x,y),t) = VectorValue(3/2*U∞*(1.0-(y/(2H))^2),0.0)
  uin(t::Real) = x -> uin(x,t)
  utop((x,y),t) = VectorValue(0.0,0.0)
  utop(t::Real) = x -> utop(x,t)
  ϕin((x,y),t) = 1.0 * (y<H)
  ϕin(t::Real) = x -> ϕin(x,t)

  # Define the finite element spaces
  @unpack order = params
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
  reffeᵩ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵤ, conformity=:H1, dirichlet_tags=["inlet","top"])
  U = TransientTrialFESpace(V, [uin,utop])
  Q = TestFESpace(Ω,reffeₚ, conformity=:L2)
  P = TrialFESpace(Q)
  Ψ = TestFESpace(Ω,reffeᵩ, conformity=:L2, dirichlet_tags=["inlet"])
  Φ = TransientTrialFESpace(Ψ,ϕin)
  X = TransientMultiFieldFESpace([U,P,Φ])
  Y = MultiFieldFESpace([V,Q,Ψ])

  # Initial solution
  xₕ₀ = interpolate_everywhere([uin(0.0),0.0,ϕin(0.0)],X(0.0))
  filename = datadir("sims","sol0")
  writevtk(Ω,filename,cellfields=["u"=>xₕ₀[1],"p"=>xₕ₀[2],"phi"=>xₕ₀[3]])

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓfp = Measure(Γfp,degree)

  # Operators
  @unpack μ,ρw,ρs,𝒟,K = params
  res(t,(u,p,ϕ),(v,q,ψ)) = ∫( (∂t(u) + (u⋅∇(u))) ⋅ v + μ*(∇(u)⊙∇(v)) - p*(∇⋅v) + q*(∇⋅u) +
                              (∂t(ϕ) + (u⋅∇(ϕ))) ⋅ ψ + 𝒟*(∇(ϕ)⊙∇(ψ)) )dΩ -
                              # ∫( K*(jump(ϕ)*mean(ψ)) + 𝒟*(mean(∇(ϕ))⋅jump(ψ⊗nfp)))dΓfp
                              ∫( K*(jump(ϕ)*mean(ψ)) )dΓfp
  jac(t,(u,p,ϕ),(du,dp,dϕ),(v,q,ψ)) = ∫( ((du⋅∇(u)) + (u⋅∇(du))) ⋅ v + μ*(∇(du)⊙∇(v)) - dp*(∇⋅v) + q*(∇⋅du) +
                              ((du⋅∇(ϕ)) + (u⋅∇(dϕ))) ⋅ ψ + 𝒟*(∇(dϕ)⊙∇(ψ)) )dΩ -
                              # ∫( K*(jump(dϕ)*mean(ψ)) + 𝒟*(mean(∇(dϕ))⋅jump(ψ⊗nfp)))dΓfp
                              ∫( K*(jump(dϕ)*mean(ψ)) )dΓfp
  jac_t(t,(u,p,ϕ),(dut,dpt,dϕt),(v,q,ψ)) = ∫( (dut) ⋅ v + (dϕt) ⋅ ψ )dΩ
  op = TransientFEOperator(res,jac,jac_t,X,Y)

  # Solver
  @unpack Δt,tf = params
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  ode_solver = ThetaMethod(nls,Δt,0.5)

  # solution
  xₕₜ = solve(ode_solver,op,xₕ₀,0.0,tf)

  # Post-processing
  filename = datadir("sims","sol")
  createpvd(filename) do pvd
    for ((uₕ,pₕ,ϕₕ),t) in xₕₜ
      pvd[t] = createvtk(Ω,filename*"_$t",cellfields=["u"=>uₕ,"p"=>pₕ,"phi"=>ϕₕ])
    end
  end

  return nothing

end

"""
NavierStokes_ConvectionDiffusion_params

This type defines a Parameters object with the default parameters for the
  NS_CD problem.
"""
@with_kw struct NavierStokes_ConvectionDiffusion_params
  H::Float64 = 1.0 # Height of the half-channel
  L::Float64 = 2.0 # Length of the channel
  μ::Float64 = 1.0 # Viscosity
  𝒟::Float64 = 1.0 # Diffusion coefficient
  K::Float64 = 1.0 # Permeability of the membrane
  ρw::Float64 = 1.0 # Density of water
  ρs::Float64 = 1.0 # Density of salt
  ne::Int = 2 # Number of elements in each direction
  order::Int = 1 # Order of the finite elements
  U∞::Float64 = 1.0 # Inlet velocity
  Δt::Float64 = 0.1 # Time step
  tf::Float64 = 0.1 # Final time
end
end
