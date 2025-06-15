export arc, circle, annulus, partial_cylinder, cylinder, partial_tube, tube, rectangle, cube
export hole_in_square_plate


function circle(; origin=(0.0,0.0), radius::Real=1.0, codimension::Int=2)
    pts, wts, θ = arc_cpts_and_wts(origin, radius, 0.0, 2π)
    space = SplineSpace(2, θ; cperiodic=[1])

    curve = GeometricMapping(Nurbs, space; codimension=codimension)
    curve.weights .= wts[1:end-1]
    curve[1].coeffs .= pts[1:end-1,1]
    curve[2].coeffs .= pts[1:end-1,2]

    return curve
end

function arc(; origin=(0.0,0.0), radius::Real=1.0, α::Real=0.0, β::Real=π/2, codimension::Int=2)
    pts, wts, θ = arc_cpts_and_wts(origin, radius, α, β)
    space = SplineSpace(2, θ)

    curve = GeometricMapping(Nurbs, space; codimension=codimension)
    curve.weights .= wts
    curve[1].coeffs .= pts[:,1]
    curve[2].coeffs .= pts[:,2]

    return curve
end

function arc_cpts_and_wts(o, r, α, β)
    # initialize
    n  = steps(β-α)
    θ  = LinRange(α, β, n+1)
    Δθ = θ[2] - θ[1]
    R  = rotation(Δθ)

    # compute weights and points of a single element
    p = cpts_reference_arc(r, Δθ)
    w = wts_reference_arc(Δθ)

    # compute weights and points of multi-element discretization
    m = 2n + 1
    wts, pts = zeros(m), zeros(m,2)
    i₁, i₂ = 1, 3
    for k in 1:n
        wts[i₁:i₂] .= w
        pts[i₁:i₂,:] .= p * rotation(θ[k])' .+ [o[1]  o[2]] # affine transformation
        pts[i₁:i₂,:] .*= w

        i₁ = i₂
        i₂ += 2
    end

    # construct knot vector
    θ = [θ...]
    m = fill(2, length(θ)); m[1] = m[end] = 3 # multiplicities

    return pts, wts, KnotVector(θ, m)
end

function steps(θ)
    return ceil(Int, θ / (2π/3))
end

# compute the control-points and weights corresponding to a quadratic
# Nurbs arc of radius 'r' and angle 'θ'
function cpts_reference_arc(r::Real, θ::Real)
    pts = zeros(3,2)
    β = θ / 2
    γ = (π / 2) - β
    pts[1,:] = r .* [1.0, 0.0]
    pts[2,:] = (r / sin(γ)) .* [cos(β), sin(β)]
    pts[3,:] = r .* [cos(2β), sin(2β)]
    return pts
end

wts_reference_arc(θ) = [1.0, sin((π-θ)/2), 1.0]

rotation(α) = [cos(α) -sin(α);
               sin(α)  cos(α)]

function annulus(;origin=(0.0,0.0), inner_radius=1.0, outer_radius=2.0, α = 0.0, β = π/2)

    # get boundary arcs
    c₁ = arc(origin=origin, radius=outer_radius, α=α, β=β)
    c₂ = arc(origin=origin, radius=inner_radius, α=α, β=β)

    # create spline spaces
    Sᵤ = c₁.space
    Sᵥ = SplineSpace(Degree(1), Interval(inner_radius, outer_radius), 1)
    
    # create mapping
    surf = GeometricMapping(Nurbs, Sᵤ ⨷ Sᵥ; codimension=2)

    # specify weights
    surf.weights[:,1]   = c₁.weights
    surf.weights[:,2]   = c₂.weights
    
    # specify x-coordinates
    surf[1].coeffs[:,1] = c₁[1].coeffs
    surf[1].coeffs[:,2] = c₂[1].coeffs

    # specify y-coordinates
    surf[2].coeffs[:,1] = c₁[2].coeffs
    surf[2].coeffs[:,2] = c₂[2].coeffs

    return surf
end

function cylinder(; origin=(0.0,0.0,0.0), radius::Real=1.0, height=1.0)
    c = circle(origin=(origin[1], origin[2]), radius=radius)

    Sᵤ = c.space
    Sᵥ = SplineSpace(Degree(1), KnotVector([0,height], [2,2]))
    surf = GeometricMapping(Nurbs, Sᵤ ⨷ Sᵥ; codimension=3)

    surf.weights[:,1]   = surf.weights[:,2]   = c.weights
    surf[1].coeffs[:,1] = surf[1].coeffs[:,2] = c[1].coeffs
    surf[2].coeffs[:,1] = surf[2].coeffs[:,2] = c[2].coeffs
    surf[3].coeffs[:,1] .= origin[3].*c.weights
    surf[3].coeffs[:,2] .= (origin[3]+height).*c.weights

    return surf
end

function partial_cylinder(; origin=(0.0,0.0,0.0), radius::Real=1.0, height=1.0, α = 0.0, β = π/2)
    c = arc(origin=(origin[1], origin[2]), radius=radius, α=α, β=β)

    Sᵤ = c.space
    Sᵥ = SplineSpace(Degree(1), KnotVector([0,height], [2,2]))
    surf = GeometricMapping(Nurbs, Sᵤ ⨷ Sᵥ; codimension=3)

    surf.weights[:,1]   = surf.weights[:,2]   = c.weights
    surf[1].coeffs[:,1] = surf[1].coeffs[:,2] = c[1].coeffs
    surf[2].coeffs[:,1] = surf[2].coeffs[:,2] = c[2].coeffs
    surf[3].coeffs[:,1] .= origin[3].*c.weights
    surf[3].coeffs[:,2] .= (origin[3]+height).*c.weights

    return surf
end

# ToDo: combine code of the two functions below
function partial_tube(; inner_radius::Real=1.0, outer_radius::Real=1.0, height=1.0, α = 0.0, β = π/2)
    t = outer_radius - inner_radius
    @assert t > 0
    c_inner = arc(radius=inner_radius, α=α, β=β)
    c_outer = arc(radius=outer_radius, α=α, β=β)

    S₁ = SplineSpace(Degree(1), KnotVector([0,t], [2,2]))
    S₂ = c_inner.space
    S₃ = SplineSpace(Degree(1), KnotVector([0,height], [2,2]))
    surf = GeometricMapping(Nurbs, S₁ ⨷ S₂ ⨷ S₃; codimension=3)

    surf.weights[1,:,1] = surf.weights[1,:,2] = c_inner.weights
    surf.weights[2,:,1] = surf.weights[2,:,2] = c_outer.weights

    surf[1].coeffs[1,:,1] = surf[1].coeffs[1,:,2] = c_inner[1].coeffs
    surf[2].coeffs[1,:,1] = surf[2].coeffs[1,:,2] = c_inner[2].coeffs
    surf[3].coeffs[1,:,1] .= 0.0
    surf[3].coeffs[1,:,2] .= height .* c_inner.weights

    surf[1].coeffs[2,:,1] = surf[1].coeffs[2,:,2] = c_outer[1].coeffs
    surf[2].coeffs[2,:,1] = surf[2].coeffs[2,:,2] = c_outer[2].coeffs
    surf[3].coeffs[2,:,1] .= 0.0
    surf[3].coeffs[2,:,2] .= height .* c_outer.weights

    return surf
end

function tube(; inner_radius::Real=1.0, outer_radius::Real=1.0, height=1.0)
    t = outer_radius - inner_radius
    @assert t > 0
    c_inner = circle(radius=inner_radius)
    c_outer = circle(radius=outer_radius)

    S₁ = SplineSpace(Degree(1), KnotVector([0,t], [2,2]))
    S₂ = c_inner.space
    S₃ = SplineSpace(Degree(1), KnotVector([0,height], [2,2]))
    surf =GeometricMapping(Nurbs, S₁ ⨷ S₂ ⨷ S₃; codimension=3)

    surf.weights[1,:,1] = surf.weights[1,:,2] = c_inner.weights
    surf.weights[2,:,1] = surf.weights[2,:,2] = c_outer.weights

    surf[1].coeffs[1,:,1] = surf[1].coeffs[1,:,2] = c_inner[1].coeffs
    surf[2].coeffs[1,:,1] = surf[2].coeffs[1,:,2] = c_inner[2].coeffs
    surf[3].coeffs[1,:,1] .= 0.0
    surf[3].coeffs[1,:,2] .= height .* c_inner.weights

    surf[1].coeffs[2,:,1] = surf[1].coeffs[2,:,2] = c_outer[1].coeffs
    surf[2].coeffs[2,:,1] = surf[2].coeffs[2,:,2] = c_outer[2].coeffs
    surf[3].coeffs[2,:,1] .= 0.0
    surf[3].coeffs[2,:,2] .= height .* c_outer.weights

    return surf
end

function rectangle(; width::Real=1.0, height::Real=1.0)
    @assert width > 0 && height > 0

    S₁ = SplineSpace(Degree(1), KnotVector([0., width], [2,2]))
    S₂ = SplineSpace(Degree(1), KnotVector([0., height], [2,2]))
    
    surf = GeometricMapping(Nurbs, S₁ ⨷ S₂; codimension=2)

    surf.weights .= 1
    
    surf[1].coeffs[1,:] .= 0.
    surf[1].coeffs[2,:] .= width

    surf[2].coeffs[:,1] .= 0.
    surf[2].coeffs[:,2] .= height
    
    return surf
end

function cube(; width::Real=1.0, height::Real=1.0, depth::Real=1.0)
    @assert depth > 0
    
    surf = rectangle(width=width, height=height)
    S₃ = SplineSpace(Degree(1), KnotVector([0., depth], [2,2]))
    
    volume = GeometricMapping(Nurbs, surf.space ⨷ S₃; codimension=3)

    volume.weights .= 1
    
    volume[1].coeffs[:,:,1] = volume[1].coeffs[:,:,2] = surf[1].coeffs
    volume[2].coeffs[:,:,1] = volume[2].coeffs[:,:,2] = surf[2].coeffs
    
    volume[3].coeffs[:,:,1] .= 0.
    volume[3].coeffs[:,:,2] .= depth
    
    return volume
end


function hole_in_square_plate()

    # create spline spaces
    Sᵤ = SplineSpace(Degree(2), Interval(0.0, pi/2), 2)
    Sᵥ = SplineSpace(Degree(2), Interval(1.0, 4.0), 1)

    # create mapping
    surf = GeometricMapping(Nurbs, Sᵤ ⨷ Sᵥ; codimension=2)

    w = (1 + 1/sqrt(2))/2
    x = sqrt(2)-1.0
    surf.weights[:] .= [ 1.0;   w;    w;   1.0;  1.0;  1.0;   1.0;  1.0;   1.0;  1.0;  1.0; 1.0]     
    surf[1].coeffs[:] .= [0.0; x; 1.0; 1.0; 0.0; 0.75; 2.5; 2.5; 0.0; 4.0; 4.0; 4.0] .* surf.weights[:]
    surf[2].coeffs[:] .= [1.0; 1.0; x; 0.0; 2.5; 2.5; 0.75; 0.0; 4.0; 4.0; 4.0; 0.0] .* surf.weights[:]

    return surf
end