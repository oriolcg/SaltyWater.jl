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
using Gridap.FESpaces: zero_free_values, interpolate!
using Gridap.Fields: meas
using Gridap.Geometry
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
    array = Bool[]
    for coordsᵢ in coords
      y_avg = 0.0
      for coordᵢ in coordsᵢ
        y_avg += coordᵢ[2]
      end
      y_avg /= length(coordsᵢ)
      push!(array, y_avg < H)
    end
    return array
  end
  coords = get_cell_coordinates(Ω)
  feed_channel_mask = is_in_feed_channel(coords)
  feed_channel_indeces = findall(feed_channel_mask)
  permeate_channel_indeces = findall(.!feed_channel_mask)
  Ωf = Interior(Ω, feed_channel_indeces)
  Ωp = Interior(Ω, permeate_channel_indeces)

  # Define the interface
  function is_in_interface(coords) # this only checks for
    array = Bool[]
    for coordsᵢ in coords
      y_avg = 0.0
      x_min = 0.0
      for coordᵢ in coordsᵢ
        y_avg += coordᵢ[2]
        x_min = max(x_min,coordᵢ[1])
      end
      y_avg /= length(coordsᵢ)
      push!(array, (y_avg ≈ H) && (x_min > 0.0))
    end
    return array
end

  # Define boundary tags
  labels_Ω = get_face_labeling(𝒯)
  add_tag_from_tags!(labels_Ω,"top",[4,6])       # assign the label "top" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[1,3,7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)

  entity_tag_left_p = num_entities(labels_Ω) + 1 # add a new tag for the interface

  # this for-loop finds all the vertices and edges in interface and assigns the new tag to them
  for d in 0:1
      face_coords = get_cell_coordinates(Grid(ReferenceFE{d}, 𝒯))
      interface  = findall(is_in_interface(face_coords))

      for i in interface
          labels_Ω.d_to_dface_to_entity[d+1][i] = entity_tag_left_p
      end
  end
  add_tag!(labels_Ω,"interface",[entity_tag_left_p])

  # Define boundaries
  Γfp = Interface(Ωf,Ωp)
  Γin = Boundary(Ω, tags="inlet")
  Γout = Boundary(Ω, tags="outlet")
  Γtop = Boundary(Ω, tags="top")
  Γb = Boundary(Ω, tags="bottom")
  nfp = get_normal_vector(Γfp)
  nb = get_normal_vector(Γb)
  nout = get_normal_vector(Γout)

  # Boundary condition
  @unpack U∞ = params
  uin((x,y),t) = VectorValue(3/2*U∞*(1.0-(y/H)^2),0.0)*(y<H) + VectorValue(0.0,0.0)*(y>=H)
  uin(t::Real) = x -> uin(x,t)
  utop((x,y),t) = VectorValue(0.0,0.0)
  utop(t::Real) = x -> utop(x,t)
  ϕin((x,y),t) = 35000 * (y<H)
  ϕin(t::Real) = x -> ϕin(x,t)
  pout((x,y)) = 6.0e3 * (y<H)

  # Define the finite element spaces
  @unpack order = params
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
  reffeᵩ = ReferenceFE(lagrangian,Float64,order-1)
  V = TestFESpace(Ω,reffeᵤ, conformity=:C0, dirichlet_tags=["inlet","top","bottom","interface"],dirichlet_masks=[(true,true),(true,true),(false,true),(true,false)])
  U = TransientTrialFESpace(V, [uin,utop,utop,utop])
  Qf = TestFESpace(Ωf,reffeₚ, conformity=:C0)
  Qp = TestFESpace(Ωp,reffeₚ, conformity=:C0)
  Pf = TrialFESpace(Qf)
  Pp = TrialFESpace(Qp)
  Ψf = TestFESpace(Ωf,reffeᵩ, conformity=:H1, dirichlet_tags=["inlet"])
  Φf = TransientTrialFESpace(Ψf,ϕin)
  Ψp = TestFESpace(Ωp,reffeᵩ, conformity=:H1, dirichlet_tags=["inlet"])
  Φp = TransientTrialFESpace(Ψp,ϕin)
  X = TransientMultiFieldFESpace([U,Pf,Pp,Φf,Φp])
  Y = MultiFieldFESpace([V,Qf,Qp,Ψf,Ψp])

  Η = TrialFESpace(V, [utop(0.0),utop(0.0),utop(0.0),utop(0.0)])
  ΗΦf = TrialFESpace(Ψf,0.0)
  ΗΦp = TrialFESpace(Ψp,0.0)

  # Initial solution
  xₕ₀ = interpolate_everywhere([uin(0.0),0.0,0.0,ϕin(0.0),ϕin(0.0)],X(0.0))
  filename = datadir("sims","sol0")
  writevtk(Ω,filename,cellfields=["u"=>xₕ₀[1],"pf"=>xₕ₀[2],"pp"=>xₕ₀[3],"phif"=>xₕ₀[4],"phip"=>xₕ₀[5]],order=order)

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΩf = Measure(Ωf,degree)
  dΩp = Measure(Ωp,degree)
  dΓfp = Measure(Γfp,degree)
  dΓb = Measure(Γb,degree)
  dΓout = Measure(Γout,degree)

  # Explicit FE functions
  global ηₙₕ = interpolate(utop(0),Η)
  global uₙₕ = interpolate(uin(0.0),U(0.0))
  global fv_u = zero_free_values(U(0.0))
  global ηϕfₙₕ = interpolate(0.0,ΗΦf)
  global ϕfₙₕ = interpolate(ϕin(0),Φf(0.0))
  global fv_ϕf = zero_free_values(Φf(0.0))
  global ηϕpₙₕ = interpolate(0.0,ΗΦp)
  global ϕpₙₕ = interpolate(ϕin(0),Φp(0.0))
  global fv_ϕp = zero_free_values(Φp(0.0))

  # Physics parameters
  @unpack μ,ρw,ρs,𝒟,K,C,T = params
  ν = μ/ρw

  # Stabilization Parameters
  c₁ = 12.0
  c₂ = 2.0
  cc = 4.0
  h2 = CellField(get_cell_measure(Ω),Ω)
  h = CellField(lazy_map(dx->dx^(1/2),get_cell_measure(Ω)),Ω)
  τₘ = 1/(c₁*ν/h2 + c₂*(meas∘uₙₕ)/h)
  τc = 0.0*cc *(h2/(c₁*τₘ))
  τₘᵩ = 1/(c₁*(𝒟)/h2)# + c₂*(meas∘uₙₕ)/h)

  # Auxiliar jump Operators
  jumpfpn(uf,up) = uf.⁺⋅nfp.⁺ - up.⁻⋅nfp.⁺
  jumpfp(uf,up) = uf.⁺ - up.⁻
  meanfp(uf,up) = uf.⁺ + up.⁻

#  mean(u))⋅(jumpfp(ϕf,ϕp)*meanfp(ψf,ψp)

  # Operators
  res(t,(u,pf,pp,ϕf,ϕp),(v,qf,qp,ψf,ψp)) = ∫( (∂t(u) + (u⋅∇(u))) ⋅ v + ν*(∇(u)⊙∇(v)) +
                              τₘ*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅u)) + τc*((∇⋅u)*(∇⋅v)) )dΩ +
                              ∫( qf*(∇⋅u) - pf*(∇⋅v) + ρw*(∂t(ϕf) + (u⋅∇(ϕf))) ⋅ ψf + 𝒟*(∇(ϕf)⊙∇(ψf)) + τₘᵩ*((∇(ϕf)'⋅u - ηϕfₙₕ)⋅(∇(ψf)'⋅u)) )dΩf +
                              ∫( qp*(∇⋅u) - pp*(∇⋅v) + ρw*(∂t(ϕp) + (u⋅∇(ϕp))) ⋅ ψp + 𝒟*(∇(ϕp)⊙∇(ψp)) + τₘᵩ*((∇(ϕp)'⋅u - ηϕpₙₕ)⋅(∇(ψp)'⋅u)) )dΩp -
                              # ∫(pf*nb⋅v)dΓb -
                              ∫( (mean(u)⋅nfp.⁺)*(jumpfp(ϕf,ϕp)*meanfp(ψf,ψp)) + 𝒟*(meanfp(∇(ϕf),∇(ϕp))⋅jumpfpn(ψf,ψp)) -
                                 1/ρw*(1/K*(mean(u)⋅nfp.⁺)*(mean(v)⋅nfp.⁺) + C*T*(jumpfp(ϕf,ϕp)*(mean(v)⋅nfp.⁺))) + ν*(jump(∇(u)⋅nfp)⋅mean(v)))dΓfp +
                              ∫( pout*nout⋅v )dΓout
  jac(t,(u,pf,pp,ϕf,ϕp),(du,dpf,dpp,dϕf,dϕp),(v,qf,qp,ψf,ψp)) = ∫( ((du⋅∇(u)) + (u⋅∇(du))) ⋅ v + μ*(∇(du)⊙∇(v)) +
                                          τₘ*((∇(u)'⋅u - ηₙₕ)⋅(∇(v)'⋅du) + (∇(du)'⋅u + ∇(u)'⋅du)⋅(∇(v)'⋅u)) +
                                          τc*((∇⋅du)*(∇⋅v)))dΩ +
                                          ∫( qf*(∇⋅du) - dpf*(∇⋅v) + ρw*((u⋅∇(dϕf)) + (du⋅∇(ϕf))) ⋅ ψf + 𝒟*(∇(dϕf)⊙∇(ψf)) +
                                          τₘᵩ*( (∇(ϕf)'⋅u - ηϕfₙₕ)⋅(∇(ψf)'⋅du) + (∇(dϕf)'⋅u + ∇(ϕf)'⋅du)⋅(∇(ψf)'⋅u)) )dΩf +
                                          ∫( qp*(∇⋅du) - dpp*(∇⋅v) + ρw*((u⋅∇(dϕp)) + (du⋅∇(ϕp))) ⋅ ψp + 𝒟*(∇(dϕp)⊙∇(ψp)) +
                                          τₘᵩ*( (∇(ϕp)'⋅u - ηϕpₙₕ)⋅(∇(ψp)'⋅du) + (∇(dϕp)'⋅u + ∇(ϕp)'⋅du)⋅(∇(ψp)'⋅u)))dΩp -
                                          # ∫( dpf*nb⋅v)dΓb -
                                          ∫( ((mean(du)⋅nfp.⁺)*(jumpfp(ϕf,ϕp)*meanfp(ψf,ψp)) + (mean(u)⋅nfp.⁺)*(jumpfp(dϕf,dϕp)*meanfp(ψf,ψp))) +
                                             𝒟*(meanfp(∇(dϕf),∇(dϕp))⋅jumpfpn(ψf,ψp)) -
                                             1/ρw*(1/K*(mean(du)⋅nfp.⁺)*(mean(v)⋅nfp.⁺) + C*T*(jumpfp(dϕf,dϕp)*(mean(v)⋅nfp.⁺))) + ν*(jump(∇(du)⋅nfp)⋅mean(v)))dΓfp
                                          # ∫( ρw*((mean(du))⋅(jumpfpn(ϕf,ϕp) + mean(u))⋅(jumpfpn(dϕf,dϕp))*meanfp(ψf,ψp)) + ρw*𝒟*(meanfp(∇(dϕf),∇(dϕp))⋅jumpfpn(ψf,ψp)) -
                                          # 1/K*((mean(du) + C*T*(jumpfpn(dϕf,dϕp)) ⋅mean(v))) )dΓfp
                                      # ∫( ρw*(mean(du)⋅nfp.⁺)*(jump(ϕ)*mean(ψ)) + ρw*(mean(u)⋅nfp.⁺)*(jump(dϕ)*mean(ψ)) +
                                      #    ρw*𝒟*(mean(∇(dϕ))⋅jump(ψ⊗nfp)) - 1/K*((mean(du) + C*T*jump(dϕ*nfp))⋅mean(v)))dΓfp
  jac_t(t,(u,pf,pp,ϕ),(dut,dpft,dppt,dϕft,dϕpt),(v,qf,qp,ψf,ψp)) = ∫( (dut) ⋅ v + ρw*(dϕft) ⋅ ψf )dΩf + ∫( (dut) ⋅ v + ρw*(dϕpt) ⋅ ψp )dΩp
  op = TransientFEOperator(res,jac,jac_t,X,Y)

  # Orthogonal projection
  aη(η,κ) = ∫( τₘ*(η⋅κ) )dΩ
  bη(κ) = ∫( τₘ*((∇(uₙₕ)'⋅uₙₕ)⋅κ) )dΩ
  aηϕf(η,κ) = ∫( (η⋅κ) )dΩf
  bηϕf(κ) = ∫( ((∇(ϕfₙₕ)'⋅uₙₕ)⋅κ) )dΩf
  aηϕp(η,κ) = ∫( (η⋅κ) )dΩp
  bηϕp(κ) = ∫( ((∇(ϕpₙₕ)'⋅uₙₕ)⋅κ) )dΩf
  op_proj = AffineFEOperator(aη,bη,Η,V)
  op_proj_ϕf = AffineFEOperator(aηϕf,bηϕf,ΗΦf,Ψf)
  op_proj_ϕp = AffineFEOperator(aηϕp,bηϕp,ΗΦp,Ψp)
  ls_proj = LUSolver()

  # Solver
  @unpack Δt,tf = params
  nls = NLSolver(show_trace=true,method=:newton,iterations=15)
  ode_solver = ThetaMethod(nls,Δt,1.0)

  # solution
  xₕₜ = solve(ode_solver,op,xₕ₀,0.0,tf)

  # Post-processing
  filename = datadir("sims","sol")
  createpvd(filename) do pvd
    for ((uₕ,pfₕ,ppₕ,ϕfₕ,ϕpₕ),t) in xₕₜ
      pvd[t] = createvtk(Ω,filename*"_$t",cellfields=["u"=>uₕ,"pf"=>pfₕ,"pp"=>ppₕ,"phif"=>ϕfₕ,"phip"=>ϕpₕ,"eta_u"=>ηₙₕ,"eta_phif"=>ηϕfₙₕ,"eta_phip"=>ηϕpₙₕ,"un"=>uₙₕ,"phifn"=>ϕfₙₕ,"phipn"=>ϕpₙₕ],order=order)
      uₙₕ = interpolate!(uₕ,fv_u,U(t))
      ϕfₙₕ = interpolate!(ϕfₕ,fv_ϕf,Φf(t))
      ϕpₙₕ = interpolate!(ϕpₕ,fv_ϕp,Φp(t))
      ηₙₕ = solve(ls_proj,op_proj)
      ηϕfₙₕ = solve(ls_proj,op_proj_ϕf)
      ηϕpₙₕ = solve(ls_proj,op_proj_ϕp)
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
  order::Int = 2 # Order of the finite elements
  U∞::Float64 = 0.06 # Inlet velocity
  Δt::Float64 = 0.1 # Time step
  tf::Float64 = 0.1 # Final time
end
end
