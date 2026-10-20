using IMAS
using Test

@testset "geometry" begin
    @testset "area of open and closed polygons" begin
        # 0.5 x 0.4 rectangle far from the origin: an open polygon used to miss the closing segment
        r = [7.0, 7.5, 7.5, 7.0]
        z = [6.0, 6.0, 6.4, 6.4]
        @test IMAS.area(r, z) ≈ 0.2
        @test IMAS.area(vcat(r, r[1]), vcat(z, z[1])) ≈ 0.2
        @test IMAS.area(reverse(r), reverse(z)) ≈ 0.2
        @test IMAS.revolution_volume(r, z) ≈ 0.2 * 2π * 7.25

        θ = range(0, 2π, 1001)
        @test IMAS.area(cos.(θ), sin.(θ)) ≈ π rtol = 1e-4                  # closed
        @test IMAS.area(cos.(θ[1:end-1]), sin.(θ[1:end-1])) ≈ π rtol = 1e-4  # open
    end

    @testset "pf_active coil area and volume" begin
        dd = IMAS.dd()
        coil = resize!(dd.pf_active.coil, 1)[1]
        element = resize!(coil.element, 1)[1]
        element.geometry.geometry_type = 2 # rectangle
        element.geometry.rectangle.r = 7.25
        element.geometry.rectangle.z = 6.2
        element.geometry.rectangle.width = 0.5
        element.geometry.rectangle.height = 0.4
        @test IMAS.area(coil) ≈ 0.2
        @test IMAS.volume(coil) ≈ 0.2 * 2π * 7.25
    end
end
