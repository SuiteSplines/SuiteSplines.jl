using IgaBase, AbstractMappings, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using LinearAlgebra, SparseArrays, StaticArrays
using TensorProductBsplines, IgaFormation, Algoim

using ImmersedSplines, Elasticity

function hole_in_plate_3d(ndofs, degree)

    # define domain and mapping
    dom = Interval(0.0,4.0) ⨱ Interval(0.0,4.0) ⨱ Interval(0.0,0.25)
    mapping = GeometricMapping(dom, (x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)
    material = Material(dim=3, youngs_modulus=1e5, poisson_ratio=0.3)

    # get analytical solutions
    benchmark = benchmark_hole_in_plate_3d(material=material);

    # partitioning into elements
    partition = CartesianProduct((d,n) -> IncreasingRange(d,n), domain(mapping), (ndofs,ndofs,3));

    # spline discretization
    trialspace_x = TensorProduct(
                    SplineSpace(degree, partition.data[1]; cleft=1:1),       # clamped x-displacement
                    SplineSpace(degree, partition.data[2]),
                    SplineSpace(degree, partition.data[3])
                );
    trialspace_y = TensorProduct(
                    SplineSpace(degree, partition.data[1]),
                    SplineSpace(degree, partition.data[2]; cleft=1:1),       # clamped y-displacement
                    SplineSpace(degree, partition.data[3])
                );
    trialspace_z = TensorProduct(
                    SplineSpace(degree, partition.data[1]),
                    SplineSpace(degree, partition.data[2]),
                    SplineSpace(degree, partition.data[3]; cleft=1:1)        # clamped z-displacemenent
                );
    displacement = Field(TensorProductBspline, (trialspace_x, trialspace_y, trialspace_z); codimension=(1,3));

    solve!(displacement=displacement, 
            mapping=mapping, 
            distancefun=benchmark.distancefun, 
            material=material, 
            traction=benchmark.stress)

    # get the stresses
    stress = CauchyStress(mapping, displacement, material);
   
    # get analytical displacement and strain. These still need to be evaluated as scalar functions.
    u = Field{1,3}( (x,y,z)->benchmark.displacement(x,y,z)[1],
                    (x,y,z)->benchmark.displacement(x,y,z)[2], 
                    (x,y,z)->benchmark.displacement(x,y,z)[3])

    σ = Field{3,3}(
            (x,y,z)->benchmark.stress(x,y,z)[1,1], 
            (x,y,z)->benchmark.stress(x,y,z)[2,1], 
            (x,y,z)->benchmark.stress(x,y,z)[3,1], 
            (x,y,z)->benchmark.stress(x,y,z)[1,2], 
            (x,y,z)->benchmark.stress(x,y,z)[2,2], 
            (x,y,z)->benchmark.stress(x,y,z)[3,2], 
            (x,y,z)->benchmark.stress(x,y,z)[1,3], 
            (x,y,z)->benchmark.stress(x,y,z)[2,3], 
            (x,y,z)->benchmark.stress(x,y,z)[3,3])

    # compue L2 error in displacement and stress
    e_displacement = l2_error(benchmark.distancefun, displacement, to=u ∘ mapping, relative=true)
    e_stress = l2_error(benchmark.distancefun, stress, to=σ ∘ mapping, relative=true)

    return e_displacement, e_stress
end

N = [10, 20, 40, 80]

E_d = zeros(length(N), 3)
E_s = zeros(length(N), 9)

count = 1
for p in 2
    for n in N

        e_displacement, e_stress = hole_in_plate_3d(n, p)
        E_d[count, :] = e_displacement
        E_s[count, :] = e_stress 

        count += 1
        println("Case degree $p and mesh parameter $n.")
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
