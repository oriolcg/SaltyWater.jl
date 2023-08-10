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
  @unpack H,L,nex,ney = params
  𝒯 = CartesianDiscreteModel((0,L,0,2H), (nex,ney))
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
  nb = get_normal_vector(Γb)

  # Boundary condition
  @unpack U∞ = params
  uin((x,y),t) = VectorValue(3/2*U∞*(1.0-(y/H)^2),0.0)*(y<H) + VectorValue(0.0,0.0)*(y>=H)
  uin(t::Real) = x -> uin(x,t)
  utop((x,y),t) = VectorValue(0.0,0.0)
  utop(t::Real) = x -> utop(x,t)
  ϕin((x,y),t) = 35000 * (y<H)
  ϕin(t::Real) = x -> ϕin(x,t)

  # Define the finite element spaces
  @unpack order = params
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
  reffeᵩ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵤ, conformity=:H1, dirichlet_tags=["inlet","top","bottom"],dirichlet_masks=[(true,true),(true,true),(false,true)])
  U = TransientTrialFESpace(V, [uin,utop,utop])
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
  dΓb = Measure(Γb,degree)

  # Operators
  @unpack μ,ρw,ρs,𝒟,K,C,T = params
  ν = μ/ρw
  res(t,(u,p,ϕ),(v,q,ψ)) = ∫( (∂t(u) + (u⋅∇(u))) ⋅ v + ν*(∇(u)⊙∇(v)) - p*(∇⋅v) + q*(∇⋅u) +
                              ρw*(∂t(ϕ) + (u⋅∇(ϕ))) ⋅ ψ + ρw*𝒟*(∇(ϕ)⊙∇(ψ)) )dΩ +
                              ∫(p*nb⋅v)dΓb -
                              ∫( ρw*(mean(u)⋅nfp.⁺)*(jump(ϕ)*mean(ψ)) + ρw*𝒟*(mean(∇(ϕ))⋅jump(ψ⊗nfp)) -
                                 1/K*((mean(u) + C*T*jump(ϕ*nfp))⋅mean(v)))dΓfp
  jac(t,(u,p,ϕ),(du,dp,dϕ),(v,q,ψ)) = ∫( ((du⋅∇(u)) + (u⋅∇(du))) ⋅ v + μ*(∇(du)⊙∇(v)) - dp*(∇⋅v) + q*(∇⋅du) +
                                          ρw*((du⋅∇(ϕ)) + (u⋅∇(dϕ))) ⋅ ψ + 𝒟*(∇(dϕ)⊙∇(ψ)) )dΩ +
                                          ∫(dp*nb⋅v)dΓb -
                                      ∫( ρw*(mean(du)⋅nfp.⁺)*(jump(ϕ)*mean(ψ)) + ρw*(mean(u)⋅nfp.⁺)*(jump(dϕ)*mean(ψ)) +
                                         ρw*𝒟*(mean(∇(dϕ))⋅jump(ψ⊗nfp)) - 1/K*((mean(du) + C*T*jump(dϕ*nfp))⋅mean(v)))dΓfp
  jac_t(t,(u,p,ϕ),(dut,dpt,dϕt),(v,q,ψ)) = ∫( (dut) ⋅ v + ρw*(dϕt) ⋅ ψ )dΩ
  op = TransientFEOperator(res,jac,jac_t,X,Y)

  # Solver
  @unpack Δt,tf = params
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  ode_solver = ThetaMethod(nls,Δt,1.0)

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
  H::Float64 = 4.0e-1 # Height of the half-channel
  L::Float64 = 1.0 # Length of the channel
  μ::Float64 = 1.0e-0 # Dynamic viscosity
  𝒟::Float64 = 1.69e-0 # Diffusion coefficient
  K::Float64 = 3.6e-0 # Permeability of the membrane
  C::Float64 = 0.2641 # Concentration vs osmotic pressure coefficient
  T::Float64 = 298.0 # Temperature
  ρw::Float64 = 1.0e0 # Density of water
  ρs::Float64 = 1.0e0 # Density of salt
  nex::Int = 2 # Number of elements in x direction
  ney ::Int = 2 # Number of elements in y direction
  order::Int = 1 # Order of the finite elements
  U∞::Float64 = 0.06 # Inlet velocity
  Δt::Float64 = 0.1 # Time step
  tf::Float64 = 0.1 # Final time
end
end
