module TestImmersedQuadrature

using Test

using IgaBase, SortedSequences, UnivariateSplines, CartesianProducts, KroneckerProducts
using IgaFormation, Algoim, ImmersedSplines
using LinearAlgebra, SparseArrays, StaticArrays

@testset "Area and perimeter quadrature" begin
    
    # circle
    phi = AlgoimCallLevelSetFunction(
            (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
            (x) -> [2.0 * x[1], 2.0 * x[2]]
            )

    # spline discretization - quarter circle
    partition = IncreasingRange(0.0,1.1,4) ⨱ IncreasingRange(0.0,1.1,5)
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]);
    
    # define regular quadrature rule
    Q = TensorProduct((d, u) -> PatchRule(d; npoints=3, method=Legendre), partition, U);
    acc_standard = ElementAccessor(testspace=U, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=false);

    # define cutcell quadrature rule
    Q = CutcellQuadratureRule(partition=partition, mapping=phi, npoints=12);
    acc_immersed = ElementAccessor(testspace=U, trialspace=U, quadrule=Q);

    # compute the volume of the cylinder
    m = prod(u -> u.p+1, U)

    # surface integrals
    Area = zeros(m,m)
    area = 0.0
    for e in Elements(partition)

        # get the element test
        etest = InsideOutsideTest(phi, e)
        
        if !is_outside(etest)
    
            # get test and trial indices
            A = TestIndices(acc_standard, e)
            B = TrialIndices(acc_standard, e)

            # case 1: element is inside
            if is_inside(etest)
                
                w = QuadratureWeights(acc_standard, e)
                u = TrialFunctions(acc_standard, e, ders=(0,0))
                v = TestFunctions(acc_standard, e, ders=(0,0))

                for k in eachindex(w)
                    Area += (v[:,k] * w[k]) * u[:,k]'
                end
                area += sum(w)

            # case 2: on the interface
            elseif on_interface(etest)
                w = QuadratureWeights(acc_immersed, e; phase=-1)
                for k in eachindex(w)
                    u = TrialFunctions(acc_immersed, e, k; ders=(0,0))
                    v = TestFunctions(acc_immersed, e, k; ders=(0,0))
                    Area += (v * w[k]) * u'
                end
                area += sum(w)
            end

        end # !is_outside(etest)
    end
    @test abs(area - π/4) < 1e-8
    @test abs(sum(Area) - π/4) < 1e-8
    
    # compute surface area of the cylinder
    Len = zeros(m, m)
    len = 0.0
    # loop over body elements (which includes the curved sides)
    for e in Elements(partition)

        # get the element test
        etest = InsideOutsideTest(phi, e)
    
        if on_interface(etest)
            w = QuadratureWeights(acc_immersed, e; phase=0)
            for k in eachindex(w)
                u = TrialFunctions(acc_immersed, e, k; ders=(0,0))
                v = TestFunctions(acc_immersed, e, k; ders=(0,0))
                Len += (v * w[k]) * u'
            end
            len += sum(w)
        end
    end
    @test abs(len - 2π/4) < 1e-8
    @test abs(sum(Len) - 2π/4) < 1e-8

end

@testset "Volume and surface boundary quadrature" begin
    
    # cylinder
    phi = AlgoimCallLevelSetFunction(
        (x) -> x[1]*x[1] + x[2]*x[2] - 1.0, 
        (x) -> [-2.0 * x[1], -2.0 * x[2], 0.0]
        )

    # spline discretization
    partition = IncreasingRange(-1.1,1.1,4) ⨱ IncreasingRange(-1.1,1.1,5) ⨱ IncreasingRange(0.0,1.0,6)
    U = SplineSpace(2, partition.data[1]) ⨷ SplineSpace(3, partition.data[2]) ⨷ SplineSpace(4, partition.data[3]);
    
    # define regular quadrature rule
    Q = TensorProduct((d, u) -> PatchRule(d; npoints=3, method=Legendre), partition, U);
    acc_standard = ElementAccessor(testspace=U, trialspace=U, quadrule=Q, incorporate_weights_in_testfuns=false);

    # define cutcell quadrature rule
    Q = CutcellQuadratureRule(partition=partition, mapping=phi, npoints=12);
    acc_immersed = ElementAccessor(testspace=U, trialspace=U, quadrule=Q);

    # compute the volume of the cylinder
    m = prod(u -> u.p+1, U)
    
    volume = 0.0
    Volume = zeros(m, m)
    for e in Elements(partition)

        # get the element test
        etest = InsideOutsideTest(phi, e)
        
        if !is_outside(etest)
    
            # get test and trial indices
            A = TestIndices(acc_standard, e)
            B = TrialIndices(acc_standard, e)

            # case 1: element is inside
            if is_inside(etest)
                
                w = QuadratureWeights(acc_standard, e)
                u = TrialFunctions(acc_standard, e, ders=(0,0,0))
                v = TestFunctions(acc_standard, e, ders=(0,0,0))

                for k in eachindex(w)
                    Volume += (v[:,k] * w[k]) * u[:,k]'
                end
                volume += sum(w)

            # case 2: on the interface
            elseif on_interface(etest)
                w = QuadratureWeights(acc_immersed, e; phase=-1)
                for k in eachindex(w)
                    u = TrialFunctions(acc_immersed, e, k; ders=(0,0,0))
                    v = TestFunctions(acc_immersed, e, k; ders=(0,0,0))
                    Volume += (v * w[k]) * u'
                end
                volume += sum(w)
            end

        end # !is_outside(etest)
    end
    @test abs(volume - π) < 1e-8
    @test abs(sum(Volume) - π) < 1e-8
    
    # compute surface area of the cylinder
    Area = zeros(m, m)
    area = 0.0
    # loop over body elements (which includes the curved sides)
    for e in Elements(partition)

        # get the element test
        etest = InsideOutsideTest(phi, e)
        
        # case 1: element is inside
        if on_interface(etest)
            w = QuadratureWeights(acc_immersed, e; phase=0)
            for k in eachindex(w)
                u = TrialFunctions(acc_immersed, e, k; ders=(0,0,0))
                v = TestFunctions(acc_immersed, e, k; ders=(0,0,0))
                Area += (v * w[k]) * u'
            end
            area += sum(w)
        end
    end
    @test abs(area - 2π) < 1e-8
    @test abs(sum(Area) - 2π) < 1e-8


    # loop over all six sides
    Area = zeros(m, m)
    area = 0.0
    for side in 1:6
        # loop over the boundary elements
        for e in Elements(restriction(partition, side))
        
            # get the element test
            etest = InsideOutsideTest(phi, e)
            
            # case 1: element is inside
            if is_inside(etest)
                w = QuadratureWeights(acc_standard, e)
                u = TrialFunctions(acc_standard, e, ders=(0,0,0))
                v = TestFunctions(acc_standard, e, ders=(0,0,0))

                for k in eachindex(w)
                    Area += (v[:,k] * w[k]) * u[:,k]'
                end
                area += sum(w)

            # case 2: element is on the interface
            elseif on_interface(etest)
                w = QuadratureWeights(acc_immersed, e; phase=-1)
                for k in eachindex(w)
                    u = TrialFunctions(acc_immersed, e, k; ders=(0,0,0))
                    v = TestFunctions(acc_immersed, e, k; ders=(0,0,0))
                    Area += (v * w[k]) * u'
                end
                area += sum(w)
            end
        end
    end
    @test abs(area - 2π) < 1e-10
    @test abs(sum(Area) - 2π) < 1e-10

end

end # module TestImmersedQuadrature