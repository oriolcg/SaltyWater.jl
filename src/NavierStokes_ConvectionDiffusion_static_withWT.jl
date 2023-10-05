module NavierStokes_ConvectionDiffusion_static_withWT
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
using Gridap.FESpaces: zero_free_values, interpolate!
using Gridap.Fields: meas
using Gridap.Geometry
using Parameters
using DrWatson

export NavierStokes_ConvectionDiffusion_static_withWT_params, solve_NSCD_static_withWT

"""
    solve_NSCD_static_withWT(params)

    This function solves the steady-state Navier-Stokes equations with convection-diffusion on
    a 2D rectangular domain using the finite element method. We assume a channel
    flow with a parabolic inflow profile and an osmosis membrane in the middle of
    the channel. The membrane is modeled as a Robin boundary condition.

    # Arguments
    - `params::Parameters`: A Parameters object containing the parameters of the problem.

    # Output
    - vtk file with the solution, stored in the `data/sims` folder.
"""
function solve_NSCD_static_withWT(params)

  # Define the domain
  @unpack H,L,nex,ney = params
  γ=2.5
  function coord_map(x)
    y = x[2]
    if y<=H/2
      yᵢ = H/2-H/2*tanh(γ*(abs(y-H/2)/(H/2)))/tanh(γ)
    else
      yᵢ = H/2+H/2*tanh(γ*(abs(y-H/2)/(H/2)))/tanh(γ)
    end
    return VectorValue(x[1],yᵢ)
  end
  𝒯 = CartesianDiscreteModel((0,L,0,H), (nex,ney),map=coord_map)
  Ω = Interior(𝒯)

  # Define boundary tags
  labels_Ω = get_face_labeling(𝒯)
  add_tag_from_tags!(labels_Ω,"membrane",[1,2,3,4,5,6])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)

  # Define boundaries
  Γout = Boundary(Ω, tags="outlet")
  Γin = Boundary(Ω, tags="inlet")
  Γₘ = Boundary(Ω, tags="membrane")
  nΓₘ = get_normal_vector(Γₘ)
  nout = get_normal_vector(Γout)
  nΓin = get_normal_vector(Γin)

  # Boundary condition
  @unpack ϕ∞,pₒ = params
  utop((x,y)) = VectorValue(0.0,0.0)
  ϕin((x,y)) = ϕ∞
  pout((x,y)) = pₒ

  # Define the finite element spaces
  @unpack order = params
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
  reffeᵩ = ReferenceFE(lagrangian,Float64,order-1)
  V = TestFESpace(Ω,reffeᵤ, conformity=:C0, dirichlet_tags=["membrane"],dirichlet_masks=[(true,false)]) #only for U, first is Ux, second is Uy
  U = TrialFESpace(V, [utop])
  Q = TestFESpace(Ω,reffeₚ, conformity=:C0)
  P = TrialFESpace(Q)
  Ψ = TestFESpace(Ω,reffeᵩ, conformity=:H1, dirichlet_tags=["inlet"])
  Φ = TrialFESpace(Ψ,ϕin)
  S = ConstantFESpace(𝒯)
  Θ = TransientTrialFESpace(S)
  X = TransientMultiFieldFESpace([U,P,Φ,Θ])
  Y = MultiFieldFESpace([V,Q,Ψ,S])
  X₀ = TransientMultiFieldFESpace([U,P,Φ])
  Y₀ = MultiFieldFESpace([V,Q,Ψ])


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓₘ = Measure(Γₘ,degree)
  dΓout = Measure(Γout,degree)
  dΓin = Measure(Γin,degree)

  # Physics parameters
  @unpack μ,ρw,𝒟,ΔP,I₀,κ,Jᵣ,Jₚ,ρₐ,Rᵣ,Vᵥ,Pᵢ = params
  Cₜ(θ,Uᵥ) = 0.1
  Uᵥ(t) = 1 # Wind velocity
  Γin_measured = ∑(∫(1.0)dΓin)
  U∞(θ)=Vᵥ*θ/Γin_measured
  uin(θ) = (U∞ ∘ θ)*CellField(x->VectorValue(6*x[2]/H*(1.0-(x[2]/H)),0.0),Γin)

  # Mesh related variables
  h = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω)),Ω)
  α = 1.0

  # Stabilization Parameters
  c₁ = 4.0
  c₂ = 2.0
  h2 = CellField(get_cell_measure(Ω),Ω)
  τₘᵩ(u) = 1/(c₁*(𝒟)/h2 + c₂*((u⋅u).^(1/2))/h)

  # Operators
  res(t,(u,p,ϕ,θ),(v,q,ψ,s)) = ∫( ρw*((u⋅∇(u))⋅v) + μ*(∇(u)⊙∇(v)) + q*(∇⋅u) - p*(∇⋅v) +
                            (u⋅∇(ϕ))⋅ψ + 𝒟*(∇(ϕ)⊙∇(ψ)) +
                            τₘᵩ(u)*((∇(ϕ)'⋅u)⋅(∇(ψ)'⋅u)) )dΩ -
                         ∫( ( nΓₘ'⋅(μ*(∇(u)⋅nΓₘ - p*nΓₘ)) ) * (v⋅nΓₘ) +
                            (ϕ*(u⋅nΓₘ))*ψ )dΓₘ +
                         ∫( (u⋅nΓₘ - ((ΔP-κ*ϕ)/I₀)) * ( nΓₘ'⋅(μ*(∇(v)⋅nΓₘ - q*nΓₘ)) ) +
                            α/h * (u⋅nΓₘ - ((ΔP-κ*ϕ)/I₀)) * (v⋅nΓₘ) )dΓₘ +
                         ∫( pout*nout⋅v )dΓout +
                         ∫(( (Jᵣ+Jₚ)*∂t(θ) - 1/2*ρₐ*Rᵣ*(Uᵥ(t))^2*Cₜ(θ,Uᵥ)+Vᵥ*(Pᵢ-p) )*s/Γin_measured +
                           (u-uin(θ)) ⋅ ( μ*∇(v)⋅nΓin - q*nΓin) + α/h * (u - uin(θ)) ⋅ v )dΓin
  res₀(θ₀) = ((u,p,ϕ),(v,q,ψ)) -> ∫( ρw*((u⋅∇(u))⋅v) + μ*(∇(u)⊙∇(v)) + q*(∇⋅u) - p*(∇⋅v) +
                              (u⋅∇(ϕ))⋅ψ + 𝒟*(∇(ϕ)⊙∇(ψ)) +
                              τₘᵩ(u)*((∇(ϕ)'⋅u)⋅(∇(ψ)'⋅u)) )dΩ -
                            ∫( ( nΓₘ'⋅(μ*(∇(u)⋅nΓₘ - p*nΓₘ)) ) * (v⋅nΓₘ) +
                              (ϕ*(u⋅nΓₘ))*ψ )dΓₘ +
                            ∫( (u⋅nΓₘ - ((ΔP-κ*ϕ)/I₀)) * ( nΓₘ'⋅(μ*(∇(v)⋅nΓₘ - q*nΓₘ)) ) +
                              α/h * (u⋅nΓₘ - ((ΔP-κ*ϕ)/I₀)) * (v⋅nΓₘ) )dΓₘ +
                            ∫( pout*nout⋅v )dΓout
  op = TransientFEOperator(res,X,Y)
  op₀(θ₀) = FEOperator(res₀(θ₀),X₀,Y₀)

  # Solver
  nls = NLSolver(show_trace=true,method=:newton,iterations=10)

  # Initial solution
  @unpack U∞₀ = params
  θ₀ = U∞₀*Γin_measured/Vᵥ
  uₕ₀,pₕ₀,ϕₕ₀ = solve(op₀(θ₀))
  xₕ₀ = interpolate_everywhere([uₕ₀,pₕ₀,ϕₕ₀,θ₀],X(0))

  # Solver
  @unpack Δt,tf = params
  ode_solver = ThetaMethod(nls,Δt,1.0)

  # solution
  xₕₜ = solve(ode_solver,op,xₕ₀,0.0,tf)

  # Post-processing
  filename = datadir("sims","sol")
  createpvd(filename) do pvd
    for ((uₕ,pₕ,ϕₕ,θₕ),t) in xₕₜ
      pvd[t] = createvtk(Ω,filename*"_$t",cellfields=["u"=>uₕ,"p"=>pₕ,"phi"=>ϕₕ,"theta"=>θₕ],order=order)
    end
  end

  return nothing

end

"""
NavierStokes_ConvectionDiffusion_params

This type defines a Parameters object with the default parameters for the
  NS_CD problem.
"""
@with_kw struct NavierStokes_ConvectionDiffusion_static_withWT_params
  H::Float64 = 1.0 # Height of the half-channel
  L::Float64 = 1.0 # Length of the channel
  μ::Float64 = 1.0e-0 # Dynamic viscosity
  𝒟::Float64 = 1.69e-0 # Diffusion coefficient
  # K::Float64 = 3.6e-0 # Permeability of the membrane
  # C::Float64 = 0.2641 # Concentration vs osmotic pressure coefficient
  # T::Float64 = 298.0 # Temperature
  ρw::Float64 = 1.0e0 # Density of water
  nex::Int = 3 # Number of elements in x direction
  ney ::Int = 3 # Number of elements in y direction
  order::Int = 2 # Order of the finite elements
  U∞₀::Float64 = 0.06 # Inlet velocity
  ϕ∞::Float64 = 35000 # Initial feed concentration
  ΔP::Float64 = 4053000.0 # Pressure drop
  I₀::Float64 = 8.41e10
  κ::Float64 = 4955.144
  pₒ::Float64 = 0.0
  Jᵣ::Float64 = 1 # Inertia rotor
  Jₚ::Float64 = 1 # Inertia pump
  ρₐ::Float64 = 1 # Air density
  Rᵣ::Float64 = 1 # Rotor radius
  Vᵥ::Float64 = 1 # Volumetric displacement HPP
  Pᵢ::Float64 = 0 # Inlet pressure HPP
  Δt::Float64 = 0.1 # Time step
  tf::Float64 = 0.1 # Final time
end
end
