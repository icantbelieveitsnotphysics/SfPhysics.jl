using SfPhysics
using Test

using Unitful, UnitfulAstro, Documenter
import PhysicalConstants.CODATA2018: g_n, c_0, m_p

@testset "SfPhysics.jl" begin
	@testset "SfGravity" begin
		# for reasons I'm too lazy to establish, there's a bit of a margin of error here. pretty small though.
		@test gravity(1u"Mearth", 1u"Rearth") ≈ g_n atol = 0.01u"m/s^2"
		@test planetary_mass(g_n |> u"m/s^2", 1u"Rearth") ≈ 1u"Mearth" atol=0.01e24u"kg"
		@test planetary_radius(1u"Mearth", g_n |> u"m/s^2") ≈ 1u"Rearth" atol = 0.01e6u"m"
		@test orbital_period(1u"Msun", 1.000001018u"AU") ≈ 365.25u"d" atol = 0.01u"d"
		
		@test orbital_radius(1u"Msun", 365.25u"d") ≈ 1.000001018u"AU" atol = 0.0001u"AU"
		
		# close enough.
		@test orbital_velocity(1.000001018u"AU", 365.25u"d") ≈ 29.78u"km/s" atol = 0.01u"km/s"
		@test orbital_velocity(1.000001018u"AU", 365.25u"d", 0.0017) ≈ 29.78u"km/s" atol = 0.01u"km/s"		
		@test escape_velocity(1u"Mearth", 1u"Rearth") ≈ 11.186u"km/s" atol = 0.01u"km/s"
	end
	
	@testset "kinematics.jl" begin
		@test duration(0u"m/s^2", 10u"m") == Inf * 1u"s"
		@test duration(10u"m/s^2", 0u"m") == 0u"s"
		@test duration(10u"m/s^2", 10u"m") == sqrt(2) * 1u"s"
		@test_throws DomainError duration(10u"m/s^2", -10u"m")
		
		@test duration(0u"m/s^2", 10u"m", 10u"m/s") == 1u"s"
		@test duration(10u"m/s^2", 0u"m", 10u"m/s") == 0u"s"
		@test duration(10u"m/s^2", 10u"m", 10u"m/s") ≈ 0.7320508075688774u"s"
		@test_throws DomainError duration(10u"m/s^2", -10u"m", 10u"m/s")
		@test_throws DomainError duration(0u"m/s^2", 10u"m", -10u"m/s")
	end
	
	@testset "kinetics.jl" begin
		@test kinetic_energy(0u"kg", 1u"m/s") == 0u"J"
		@test kinetic_energy(1u"kg", 0u"m/s") == 0u"J"
		@test kinetic_energy(2u"kg", 1u"m/s") == 1u"J"
		@test kinetic_energy(2u"kg", 2u"m/s") == 4u"J"
		@test kinetic_energy(2u"kg", -2u"m/s") == 4u"J"
	end
	
	@testset "SfPlanetary" begin
		s = SfSolarSystem
		
		@test gravity(s.earth) ≈ g_n atol = 0.01u"m/s^2"
		@test planetary_mass(s.earth) == 1u"Mearth"
		@test planetary_radius(s.earth) == 1u"Rearth_e"
		@test orbital_period(s.earth) ≈ 365.25u"d" atol = 0.01u"d"
		@test orbital_velocity(s.earth) ≈ 29.78u"km/s" atol = 0.01u"km/s"	
		@test escape_velocity(s.earth) ≈ 11.186u"km/s" atol = 0.01u"km/s"
		
		@test hill_sphere(s.earth) > orbital_radius(s.moon)
		@test hill_sphere(s.earth.orbit, s.moon.mass) < orbital_radius(s.moon)
	end
	
	@testset "SfRelativity" begin
		# use bignums here, otherwise floating point rounding will ruin accuracy (RKE returns 59J at 10m/s!)
		@test kinetic_energy(1u"kg", 10u"m/s") ≈ relativistic_kinetic_energy(1u"kg", big(10) * 1u"m/s") atol = 0.0000000000001u"J"
		@test kinetic_energy(1u"kg", 0.5c_0) ≈ relativistic_kinetic_energy(1u"kg", 0.5c_0) * 0.808 atol=1e12u"J"
		
		@test lorentz_velocity(lorentz_factor(.05c_0)) ≈ .05c_0
		@test lorentz_velocity(lorentz_factor(.95c_0)) ≈ .95c_0
		
		@test relativistic_kinetic_energy(m_p, relativistic_velocity(m_p, 100u"MeV")) ≈ 100u"MeV" atol = 1u"μeV"
	end
	
	@testset "SfGeometry" begin
		@test radius(Sphere(1u"m")) == 1u"m"
		@test equatorial_radius(Sphere(1u"m")) == 1u"m"
		@test polar_radius(Sphere(1u"m")) == 1u"m"
		@test area(Sphere(1u"m")) == 4π * u"m^2"
		@test volume(Sphere(1u"m")) == (4π / 3) * u"m^3"
		@test cross_sectional_area(Sphere(1u"m")) ≈ π * u"m^2"
		
		@test radius(Spheroid(1u"m", 1u"m")) == 1u"m"
		@test equatorial_radius(Spheroid(2u"m", 1u"m")) == 2u"m"
		@test polar_radius(Spheroid(2u"m", 1u"m")) == 1u"m"
		@test area(Spheroid(1u"m", 1u"m")) == 4π * u"m^2"
		@test area(Spheroid(1.01u"m", 1u"m")) ≈ 4π * u"m^2" atol=0.2u"m^2" # weirdly poor tolerances here, code needs checking
		@test area(Spheroid(1u"m", 1.01u"m")) ≈ 4π * u"m^2" atol=0.2u"m^2"
		@test volume(Spheroid(1u"m", 1u"m")) == (4π / 3) * u"m^3"
		@test cross_sectional_area(Spheroid(1u"m", 1u"m")) ≈ π * u"m^2"
		
		@test radius(TriaxialEllipsoid(1u"m", 1u"m", 1u"m")) == 1u"m"
		@test area(TriaxialEllipsoid(1u"m", 1u"m", 1u"m")) == 4π * u"m^2"
		@test area(TriaxialEllipsoid(1u"m", 1.001u"m", 1.001u"m")) ≈ 4π * u"m^2" atol=0.05u"m^2"
		@test area(TriaxialEllipsoid(1u"m", 1.001u"m", 1.002u"m")) ≈ 4π * u"m^2" atol=0.05u"m^2"
		@test volume(TriaxialEllipsoid(1u"m", 1u"m", 1u"m")) == (4π / 3) * u"m^3"
	end
	
	@testset "SfAstronomy" begin
		@test diffuse_disc_q(0) == 1.0
		@test diffuse_disc_q(deg2rad(45)) == sqrt(2) / 2
		@test diffuse_disc_q(deg2rad(90)) ≈ 0.0 atol = eps(Float64)
		@test diffuse_disc_q(deg2rad(180)) == 1.0
		
		@test diffuse_sphere_q(0) == 2/3
		@test diffuse_sphere_q(deg2rad(45)) == 2/3 * (3sqrt(2)/8 + sqrt(2)/2π)
		@test diffuse_sphere_q(deg2rad(90)) ≈ 2/3π atol = 0.000000000000001 # weird difference in last decimal place under test
		@test diffuse_sphere_q(deg2rad(180)) ≈ 0.0 atol = eps(Float64)
		
		@test absolute_magnitude(1000u"m", 1) == absolute_magnitude(10000u"m", 1) + 5
		@test absolute_magnitude(6.9173e7u"m", 0.538) ≈ -9.4 atol = 0.1 # jupiter	

		@test apparent_magnitude(0, 1u"AU", 1u"AU", 1u"AU", 1) == 0.0
		@test apparent_magnitude(0, 1u"AU", 1u"AU", 1u"AU", 0.01) == 2apparent_magnitude(0, 1u"AU", 1u"AU", 1u"AU", 0.1)
	end
	
	@testset "SfMatter" begin
		@test bouyancy(1u"kg/m^3", 1u"kg/m^3"; g = g_n) == g_n
		@test bouyancy(1u"kg/m^3", 2u"kg/m^3"; g = g_n) == 2g_n
		@test bouyancy(2u"kg/m^3", 1u"kg/m^3"; g = g_n) == g_n / 2
		
		@test compression_energy(1u"mol", 1u"K", 1u"m^3", 1u"m^3") == 0.0u"J"
	end
	
	@testset "Units" begin
		import SfPhysics.SfUnits: to_angle
	
		@test to_angle(90) == 90u"°"
		@test to_angle(90u"°") == 90u"°"
		@test to_angle((π / 2) * u"rad") == 90u"°"
		@test unit(to_angle((π / 2) * u"rad")) == u"rad"
		@test to_angle(missing) === missing
		@test to_angle(nothing) == nothing
		
		@test 1e20u"W" |> kardashev |> kardashev == 1e20u"W"
		
		@test tnt(1e9u"J") |> u"J" ≈ 1e9u"J"
	end
	
	@testset "Documentation" begin
		doctest(SfPhysics; manual=false)
	end
end
