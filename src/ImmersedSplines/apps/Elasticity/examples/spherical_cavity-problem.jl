using Revise

using IgaBase, AbstractMappings, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using LinearAlgebra, SparseArrays, StaticArrays
using TensorProductBsplines, IgaFormation, Algoim

using ImmersedSplines, Elasticity

# model initialization
dom = Interval(0.0,4.0) ⨱ Interval(0.0,4.0) ⨱ Interval(0.0,4.0)
mapping = GeometricMapping(dom, (x,y,z) -> x, (x,y,z) -> y, (x,y,z) -> z)
material = Material(dim=3, youngs_modulus=1e5, poisson_ratio=0.3)

# get analytical solutions
benchmark = benchmark_spherical_cavity(material=material);

# partitioning into elements
partition = CartesianProduct((d,n) -> IncreasingRange(d,n), domain(mapping), (16,16,16));

# spline discretization
degree = 2
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
            
using GLMakie
GLMakie.activate!()

fig = Figure(resolution = (1200, 800))
ax = LScene(fig[1, 1], show_axis = true)

# get visualization points in parametric space
X = CartesianProduct((u,m) -> IncreasingRange(u, m), dom, (201, 201, 201));

r = (-0.5, 20.5) # σ_zz
for dir in 1:3
    for comp in 1:2

        side = comp + (dir-1)*2

        # get the evalation points at this boundary
        x = restrict_to(X, side=side)

        # evaluate coordinates, displacement, and stress
        @evaluate Y = mapping(x);
        @evaluate S = stress(x);

        # plot stress
        x, y, z = map(squeeze, Y.data);
        s = squeeze(S.data[3,3]);
        if comp==2 # invert normal such that Makie does shading on correct side
            x = x'; y = y'; z = z'; s = s'
        end
        GLMakie.surface!(ax, x, y, -z; color = s, colormap = :jet, colorrange=r)        
        GLMakie.surface!(ax, -x', y', -z'; color = s', colormap = :jet, colorrange=r)
        GLMakie.surface!(ax, -x, -y, -z; color = s, colormap = :jet, colorrange=r)
        GLMakie.surface!(ax, x', -y', -z'; color = s', colormap = :jet, colorrange=r)

        # GLMakie.surface!(ax, x, y, z; color = s, colormap = :jet, colorrange=r)        
        GLMakie.surface!(ax, -x, y, z; color = s, colormap = :jet, colorrange=r)
        GLMakie.surface!(ax, -x', -y', z'; color = s', colormap = :jet, colorrange=r)
        GLMakie.surface!(ax, x, -y, z; color = s, colormap = :jet, colorrange=r)
    end
end

# mask the circilar hole with a disk on bottom and top
r = 1.0; θ = 0.0:pi/180:2*pi; β = 0.0:pi/180:pi/2 
x, y, z = r * sin.(β) * cos.(θ)', r * sin.(β) * sin.(θ)', r * cos.(β) * ones(size(θ))';
GLMakie.surface!(ax, x, y, z; colorrange = (-1,0), highclip=(:black, 1), shading=true, transparency=false)

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
l2_error(benchmark.distancefun, displacement, to=u ∘ mapping, relative=true)
l2_error(benchmark.distancefun, stress, to=σ ∘ mapping, relative=true)