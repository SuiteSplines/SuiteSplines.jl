function hole_in_plate_2d(ndofs, order)

    # define domain and mapping
    dom = Interval(0.0,4.0) ⨱ Interval(0.0,4.0)
    mapping = GeometricMapping(dom, (x,y) -> x, (x,y) -> y)
    material = Material(dim=2, youngs_modulus=1e5, poisson_ratio=0.3)

    # get analytical solutions
    benchmark = benchmark_hole_in_plate_2d(material=material);

    # partitioning into elements
    partition = CartesianProduct((d,n) -> IncreasingRange(d,n), domain(mapping), (ndofs,ndofs))

    # spline discretization
    trialspace_x = TensorProduct(
                    SplineSpace(order, partition.data[1]; cleft=1:1),       # clamped left part
                    SplineSpace(order, partition.data[2])
                );
    trialspace_y = TensorProduct(
                    SplineSpace(order, partition.data[1]),
                    SplineSpace(order, partition.data[2]; cleft=1:1)        # clamped bottom
                );

    displacement = Field(TensorProductBspline, (trialspace_x, trialspace_y); codimension=(1,2));

    solve!(displacement=displacement, 
            mapping=mapping, 
            distancefun=benchmark.distancefun, 
            material=material, 
            traction=benchmark.stress)

    # get the stresses
    stress = CauchyStress(mapping, displacement, material);

    # get analytical displacement and strain. These still need to be evaluated as scalar functions.
    u = Field{1,2}((x,y)->benchmark.displacement(x,y)[1], (x,y)->benchmark.displacement(x,y)[2])
    σ = Field{2,2}((x,y)->benchmark.stress(x,y)[1,1], (x,y)->benchmark.stress(x,y)[2,1], (x,y)->benchmark.stress(x,y)[1,2], (x,y)->benchmark.stress(x,y)[2,2])

    # compue L2 error in displacement and stress
    e_displacement = l2_error(benchmark.distancefun, displacement, to=u ∘ mapping, relative=true)
    e_stress = l2_error(benchmark.distancefun, stress, to=σ ∘ mapping, relative=true)

    return e_displacement, e_stress
end

N = [10, 20, 40, 80, 160]

E_d = zeros(length(N), 2)
E_s = zeros(length(N), 4)

count = 1
for p in 2
    for n in N

        e_displacement, e_stress = hole_in_plate_2d(n, p)
        E_d[count, :] = e_displacement
        E_s[count, :] = e_stress 

        count += 1
    end
end

using GLMakie
GLMakie.activate!()

fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1], xscale=log10, yscale=log10)

scatter!(N, E_d[:,2])
lines!(N, E_d[:,2])

scatter!(N, E_s[:,1])
lines!(N, E_s[:,1])
