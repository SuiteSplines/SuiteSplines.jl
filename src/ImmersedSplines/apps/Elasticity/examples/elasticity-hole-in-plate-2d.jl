using IgaBase, AbstractMappings, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using LinearAlgebra, SparseArrays, StaticArrays
using TensorProductBsplines, IgaFormation, Algoim

using ImmersedSplines, Elasticity

# define domain and mapping
dom = Interval(0.0,4.0) ⨱ Interval(0.0,4.0)
mapping = GeometricMapping(dom, (x,y) -> x, (x,y) -> y)
material = Material(dim=2, youngs_modulus=1e5, poisson_ratio=0.3)

# get analytical solutions
benchmark = benchmark_hole_in_plate_2d(material=material);

# partitioning into elements
partition = CartesianProduct((d,n) -> IncreasingRange(d,n), domain(mapping), (41,41))

# spline discretization
order = 2
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

using GLMakie
GLMakie.activate!()

fig = Figure(resolution = (1200, 800))
ax = LScene(fig[1, 1], show_axis = true)

# get visualization points in parametric space
X = visualization_grid(mapping, density=(201,201));

# evaluate coordinates, displacement, and stress
@evaluate Y = mapping(X);
@evaluate Δu = displacement(X);
@evaluate S = stress(X);

x, y = map(Matrix, Y.data);
z = zeros(size(x));
s = Matrix(S.data[1,1]);
GLMakie.surface!(ax, x, y, z; color = s, colormap = :jet, colorrange=(-10,30))
GLMakie.surface!(ax, -x', y', z'; color = s', colormap = :jet, colorrange=(-10,30))
GLMakie.surface!(ax, -x, -y, z; color = s, colormap = :jet, colorrange=(-10,30))
GLMakie.surface!(ax, x', -y', z'; color = s', colormap = :jet, colorrange=(-10,30))

# mask the circilar hole with a disk on bottom and top
r = 0.0:1.0; theta = 0.0:pi/180:2*pi
x, y = r' .* cos.(theta), r' .* sin.(theta);
z = fill(0.01, size(x))
GLMakie.surface!(ax, x, y, z; color = :white)
GLMakie.surface!(ax, x, y, -z; color = :white)

# get analytical displacement and strain. These still need to be evaluated as scalar functions.
u = Field{1,2}((x,y)->benchmark.displacement(x,y)[1], (x,y)->benchmark.displacement(x,y)[2])
σ = Field{2,2}((x,y)->benchmark.stress(x,y)[1,1], (x,y)->benchmark.stress(x,y)[2,1], (x,y)->benchmark.stress(x,y)[1,2], (x,y)->benchmark.stress(x,y)[2,2])

# compue L2 error in displacement and stress
l2_error(benchmark.distancefun, displacement, to=u ∘ mapping, relative=true)
l2_error(benchmark.distancefun, stress, to=σ ∘ mapping, relative=true)