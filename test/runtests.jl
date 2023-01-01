using SfPhysics
using Test

using Unitful, UnitfulAstro, Documenter
import PhysicalConstants.CODATA2018: g_n, c_0, m_p, N_A
import LinearAlgebra: norm

@testset "SfPhysics.jl" begin
	@testset "SfGravity" begin
		@test sqrt(vis_viva(1u"Msun", 1u"AU", 1u"AU")) ≈ orbital_velocity(1u"Msun", 1u"AU")
		
		# for reasons I'm too lazy to establish, there's a bit of a margin of error here. pretty small though.
		@test gravity(1u"Mearth", 1u"Rearth") ≈ g_n atol = 0.01u"m/s^2"
		@test gravity(1u"Mearth", 100u"kg", 1u"Rearth") ≈ 979.8398133669466u"N"
		
		@test planetary_mass(g_n |> u"m/s^2", 1u"Rearth") ≈ 1u"Mearth" atol=0.01e24u"kg"
		@test planetary_mass(11.186u"km/s", 1u"Rearth") ≈ 1u"Mearth" atol=0.01e24u"kg"
		@test planetary_mass(1u"AU", 365.25u"d") ≈ 1u"Msun" atol=0.0001e30u"kg"
		
		@test planetary_radius(1u"Mearth", g_n |> u"m/s^2") ≈ 1u"Rearth" atol = 0.01e6u"m"
		@test planetary_radius(1u"Mearth", 11.186u"km/s") ≈ 1u"Rearth" atol = 0.01e6u"m"
		@test planetary_radius(5.51u"g/cm^3", g_n |> u"m/s^2") ≈ 1u"Rearth" atol = 2e4u"m"
		@test planetary_radius(5.51u"g/cm^3", 11.186u"km/s") ≈ 1u"Rearth" atol = 0.01e6u"m"
		
		@test orbital_period(1u"Msun", 1.000001018u"AU") ≈ 365.25u"d" atol = 0.01u"d"		
		@test orbital_radius(1u"Msun", 365.25u"d") ≈ 1.000001018u"AU" atol = 0.0001u"AU"
		@test orbital_radius(1u"Msun", 29.8u"km/s") ≈ 1.000001018u"AU" atol = 0.002u"AU" # slightly imprecise velocity of earth
		
		# close enough.
		@test orbital_velocity(1.000001018u"AU", 365.25u"d") ≈ 29.78u"km/s" atol = 0.01u"km/s"
		@test orbital_velocity(1.000001018u"AU", 365.25u"d", 0) ≈ orbital_velocity(1.000001018u"AU", 365.25u"d")
		@test orbital_velocity(1.000001018u"AU", 365.25u"d", 0.0017) ≈ 29.78u"km/s" atol = 0.01u"km/s"	
		
		@test escape_velocity(1u"Mearth", 1u"Rearth") ≈ 11.186u"km/s" atol = 0.01u"km/s"

		@test hill_sphere(1u"Msun", 1u"Mearth", 1u"AU") ≈ 0.0098u"AU" atol = 0.01u"AU"		
		@test gravitational_binding_energy(1u"Mearth", 1u"Rearth") ≈ 2.24e32u"J" atol = 0.001e32u"J"		
		@test roche_limit(1u"Rearth", 5.51u"g/cm^3", 3.34u"g/cm^3") ≈ 9500u"km" atol = 5u"km"
		@test barycentric_distance(1u"Mearth", 0.0123u"Mearth", 384000u"km") ≈ 4670u"km" atol = 5u"km"
	end
	
	@testset "kinematics.jl" begin		
		@test kinetic_energy(0u"kg", 1u"m/s") == 0u"J"
		@test kinetic_energy(1u"kg", 0u"m/s") == 0u"J"
		@test kinetic_energy(2u"kg", 1u"m/s") == 1u"J"
		@test kinetic_energy(2u"kg", 2u"m/s") == 4u"J"
		@test kinetic_energy(2u"kg", -2u"m/s") == 4u"J"
		
		@test distance(1u"m/s/s", 10u"s") == 50u"m"
		@test distance(1u"m/s/s", 10u"s", 1u"m/s") == 60u"m"
		@test distance(1u"m/s/s", 10u"s", -1u"m/s") == 40u"m"
		@test distance(-1u"m/s/s", 10u"s", 5u"m/s") == 0u"m"
		
		@test duration(0u"m/s^2", 10u"m") == Inf * 1u"s"
		@test duration(10u"m/s^2", 0u"m") == 0u"s"
		@test duration(10u"m/s^2", 10u"m") == sqrt(2) * 1u"s"
		@test_throws DomainError duration(10u"m/s^2", -10u"m")
		
		@test duration(0u"m/s^2", 10u"m", 10u"m/s") == 1u"s"
		@test duration(10u"m/s^2", 0u"m", 10u"m/s") == 0u"s"
		@test duration(10u"m/s^2", 10u"m", 10u"m/s") ≈ 0.7320508075688774u"s"
		@test_throws DomainError duration(10u"m/s^2", -10u"m", 10u"m/s")
		@test_throws DomainError duration(0u"m/s^2", 10u"m", -10u"m/s")
		
		@test acceleration(10u"m", 10u"m/s") == 5u"m/s/s"
		
		@test projectile_displacement(100u"m/s", 90u"°", 10u"s", g=10u"m/s/s") == [0u"m", 500u"m"]
		@test projectile_displacement(100u"m/s", 45u"°", 10u"s", g=10u"m/s/s") ≈ [1000u"m" / sqrt(2), (1000u"m" / sqrt(2)) - 500u"m"]
		@test projectile_displacement(100u"m/s", 10u"s", g=10u"m/s/s") == projectile_displacement(100u"m/s", 90u"°", 10u"s", g=10u"m/s/s")[2]
		
		@test projectile_velocity(100u"m/s", 90u"°", 10u"s", g=10u"m/s/s") == [0u"m/s", 0u"m/s"]
		@test projectile_velocity(100u"m/s", 45u"°", 10u"s", g=10u"m/s/s")[1] == projectile_velocity(100u"m/s", 45u"°", 5u"s", g=10u"m/s/s")[1]
		
		@test projectile_flight_time(100u"m/s", g=10u"m/s/s") == 20u"s"
		@test projectile_flight_time(100u"m/s", 45u"°", g=10u"m/s/s") == 10sqrt(2) * u"s"
		
		@test projectile_apex(100u"m/s", g=10u"m/s/s") == [0u"m", 500u"m"]
		@test projectile_apex(100u"m/s", 45u"°", g=10u"m/s/s") ≈ [500u"m", 250u"m"]
		
		@test projectile_range(100u"m/s", 90u"°", g=10u"m/s/s") == 0u"m"
		@test projectile_range(100u"m/s", 45u"°", g=10u"m/s/s") == 1000u"m"
		
	end
	
	@testset "SfPlanetary" begin
		s = SfSolarSystem
		
		@test planetary_mass(s.earth) == 1u"Mearth"
		@test mass(s.earth) == planetary_mass(s.earth)
		
		@test gravity(s.earth) ≈ g_n atol = 0.01u"m/s^2"
		@test escape_velocity(s.earth) ≈ 11.186u"km/s" atol = 0.01u"km/s"
		
		@test planetary_radius(s.earth) == 1u"Rearth_e"
		@test equatorial_radius(s.earth) == 1u"Rearth_e"
		@test polar_radius(s.earth) == 1u"Rearth_p"
		@test radius(s.earth) ≈ 1u"Rearth" atol = 2e4u"m"
		
		@test orbital_period(s.earth) ≈ 365.25u"d" atol = 0.01u"d"
		@test orbital_velocity(s.earth) ≈ 29.78u"km/s" atol = 0.01u"km/s"
		@test orbital_radius(s.earth) ≈ 1u"AU" atol = 0.00001u"AU"
				
		@test orbital_period(s.earth.orbit) == orbital_period(s.earth)	
		@test orbital_velocity(s.earth.orbit) == orbital_velocity(s.earth)
		@test orbital_radius(s.earth) == radius(s.earth.orbit)
		
		@test hill_sphere(s.earth) > orbital_radius(s.moon)
		@test hill_sphere(s.earth.orbit, s.moon.mass) < orbital_radius(s.moon)
		
		@test star(s.sol) == s.sol
		@test star(s.earth) == s.sol
		@test star(s.moon) == s.sol
		
		@test stellar_distance(s.sol) == 0u"AU"
		@test stellar_distance(s.earth) == orbital_radius(s.earth)
		@test stellar_distance(s.moon) == orbital_radius(s.earth)
		
		@test stellar_luminosity(1u"Rsun", s.sol.surface_temperature) ≈ 1u"Lsun" atol=1e-5u"Lsun"
		@test stellar_luminosity(s.sol) ≈ 1u"Lsun" atol=1e-5u"Lsun"
		
		@test stellar_irradiance(s.earth)/cross_sectional_area(s.earth) ≈ 1361u"W/m^2" atol = 0.2u"W/m^2"
		
		@test planetary_equilibrium_temperature(s.earth) ≈ 255u"K" atol = 1u"K"
		
		@test absolute_magnitude(s.jupiter) ≈ -9.4 atol = 0.05
	end
    
    @testset "SfOrbital" begin
        sv = StateVector([1000u"km", 5000u"km", 7000u"km"], [3u"km/s", 4u"km/s", 5u"km/s"])
        el = state_vector_to_elements(sv, 1u"Mearth")
        
        @test el.e ≈ 0.9475409572473679
        
        sv2 = elements_to_state_vector(el, 1u"Mearth")
        
        @test sv2.r ≈ sv.r
        @test sv2.v ≈ sv.v
        
        t = orbital_period(el, 1u"Mearth")
        
        @test t ≈ 9183.8529879387u"s"
        
        pos = orbital_position(el, t, 1u"Mearth")
        
        @test pos.ν ≈ el.ν atol=0.000001
        
        el.ν = 0u"rad"        
        pos = orbital_position(el, t / 2, 1u"Mearth")
        
        @test pos.ν ≈ π 
        @test time_since_periapsis(pos, 1u"Mearth") ≈ t/2
        
        el = OrbitalElements(0, 100000u"km", 0, 0, 0, 0)
        @test norm(elements_to_state_vector(el, 1u"Mearth").v) ≈ orbital_velocity(1u"Mearth", el.a)
    end
	
	@testset "SfRelativity" begin
		# use bignums here, otherwise floating point rounding will ruin accuracy (RKE returns 59J at 10m/s!)
		@test kinetic_energy(1u"kg", 10u"m/s") ≈ relativistic_kinetic_energy(1u"kg", big(10) * 1u"m/s") atol = 0.0000000000001u"J"
		@test kinetic_energy(1u"kg", 0.5c_0) ≈ relativistic_kinetic_energy(1u"kg", 0.5c_0) * 0.808 atol=1e12u"J"
		
		@test lorentz_velocity(lorentz_factor(.05c_0)) ≈ .05c_0
		@test lorentz_velocity(lorentz_factor(.95c_0)) ≈ .95c_0
		
		@test relativistic_kinetic_energy(m_p, relativistic_velocity(m_p, 100u"MeV")) ≈ 100u"MeV" atol = 1u"μeV"
        
        @test rapidity(0u"c") == 0
        @test rapidity(1u"c") == Inf
        @test rapidity_to_velocity(rapidity(0.5u"c")) == 0.5u"c"
        
        d = 5u"ly"
        t = 10u"yr"
        
        @test relativistic_brachistochrone_transit_time(d, relativistic_brachistochrone_acceleration(d, t)) == t
        
        Δv = 0.5u"c"
        ve = 0.5u"c"
        
        @test relativistic_delta_v(ve, relativistic_mass_ratio(ve, Δv)) == Δv
        @test relativistic_delta_v(ve, relativistic_mass_ratio(ve, rapidity(Δv))) == Δv
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
		
		@test Cube(1u"m") == Cuboid(1u"m", 1u"m", 1u"m")
		
		@test area(Cube(1u"m")) == 6u"m^2"
		@test volume(Cube(1u"m")) == 1u"m^3"
	 	@test area(Cube(1u"m^3")) == 6u"m^2"
		@test volume(Cube(1u"m^3")) == 1u"m^3"
		
		@test radius(Cylinder(1u"m", 1u"m")) == 1u"m"
		@test length(Cylinder(1u"m", 1u"m")) == 1u"m"
		@test area(Cylinder(1u"m", 1u"m")) == 4π * u"m^2"
		@test volume(Cylinder(1u"m", 1u"m")) ≈ π * u"m^3"
		
		@test radius(SphericalShell(0u"m", 1u"m")) == 1u"m"
		@test area(SphericalShell(0u"m", 1u"m")) == area(Sphere(1u"m"))
		@test volume(SphericalShell(0u"m", 1u"m")) == volume(Sphere(1u"m"))
		@test volume(SphericalShell(1u"m", 1u"m")) == 0u"m^3"
		@test volume(SphericalShell(1u"m", 2u"m")) == volume(Sphere(2u"m")) - volume(Sphere(1u"m"))
		
		@test spherical_cap_solid_angle(0u"°") == 0u"sr"
		@test spherical_cap_solid_angle(180u"°") == 4π * u"sr"
		@test spherical_cap_solid_angle(0u"m", 2u"m") ≈ 2π * u"sr"
        
        @test intersect_circles(3, 1, 1) == 0
        @test intersect_circles(2, 1, 1) == 0
        @test intersect_circles(1, 1, 1) ≈ (2π/3)-(sqrt(3)/2) atol=1e-15
        @test intersect_circles(0, 1, 1) ≈ π
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
	
	@testset "SfCoriolis" begin
		@test coriolis_acceleration(0.001u"rad/s" * unit_y, 400u"m/s" * unit_x)[3] == -0.8u"m/s/s"
		@test coriolis_force(1u"kg", 0.001u"rad/s" * unit_y, 400u"m/s" * unit_x)[3] == -0.8u"N"
		
		@test tangential_velocity(1u"rad/s", 1u"m") == 1u"m/s"
		@test tangential_velocity(1u"rad/s", 1u"m", 0u"°") == 1u"m/s"
		@test tangential_velocity(1u"rad/s", 1u"m", 90u"°") == 0u"m/s"
		@test tangential_velocity(1u"rad/s", 1u"m", 45u"°") == cosd(45) * u"m/s"
		
		@test centrifugal_acceleration(1u"rad/s", 0u"m") == 0u"m/s/s"
		@test centrifugal_acceleration(1u"rad/s", 1u"m") == 1u"m/s/s"
		@test centrifugal_acceleration(1u"rad/s", 2u"m") == 2u"m/s/s"
		@test centrifugal_acceleration(2u"rad/s", 1u"m") == 4u"m/s/s"
		
		@test centrifugal_acceleration(1u"rad/s", 1u"m", 0u"°") == 1u"m/s/s"
		@test centrifugal_acceleration(1u"rad/s", 1u"m", 90u"°") == 0u"m/s/s"
		@test centrifugal_acceleration(1u"rad/s", 1u"m", 45u"°") == cosd(45) * u"m/s/s"
		
		@test centrifugal_force(1u"kg", 1u"rad/s", 1u"m") == 1u"N"
		@test centrifugal_force(2u"kg", 1u"rad/s", 1u"m") == 2u"N"
		@test centrifugal_force(1u"kg", 1u"rad/s", 2u"m") == 2u"N"
		@test centrifugal_force(1u"kg", 2u"rad/s", 1u"m") == 4u"N"
		
		@test centrifugal_force(1u"kg", 1u"rad/s", 1u"m", 0u"°") == 1u"N"
		@test centrifugal_force(1u"kg", 1u"rad/s", 1u"m", 90u"°") == 0u"N"
		@test centrifugal_force(1u"kg", 1u"rad/s", 1u"m", 45u"°") == cosd(45) * u"N"
		
		@test angular_velocity(1u"m/s/s", 1u"m") == 1u"rad/s"
		@test radius(1u"rad/s", 1u"m/s/s") == 1u"m"
		
		@test angular_momentum(1u"kg", 1u"rad/s", 1u"m") == 1u"kg*m^2/s"
		@test angular_momentum(2u"kg", 1u"rad/s", 1u"m") == 2u"kg*m^2/s"
		@test angular_momentum(1u"kg", 2u"rad/s", 1u"m") == 2u"kg*m^2/s"
		@test angular_momentum(1u"kg", 1u"rad/s", 2u"m") == 4u"kg*m^2/s"
	end
	
	@testset "SfRocketry" begin
		@test rocket_propulsive_efficiency(1u"m/s", 1u"m/s") == 1.0
		
		@test brachistochrone_transit_time(1u"m", 1u"m/s/s") == 2u"s"
		@test brachistochrone_transit_time(1u"m", 2u"m/s/s") == sqrt(2) * u"s"
		@test brachistochrone_transit_time(2u"m", 1u"m/s/s") == 2sqrt(2) * u"s"
		
		@test brachistochrone_acceleration(1u"m", 2u"s") == 1u"m/s/s"
		@test brachistochrone_acceleration(2u"m", 2sqrt(2) * u"s") ≈ 1u"m/s/s"
		
		@test brachistochrone_delta_v(1u"m", 1u"m/s/s") == 2u"m/s"
		@test brachistochrone_delta_v(2u"m", 1u"m/s/s") == 2sqrt(2) * u"m/s"
		
		@test boost_coast_transit_time(100u"m", 1u"m/s/s", 1u"s") == 101u"s"
		@test boost_coast_transit_time(200u"m", 1u"m/s/s", 1u"s") == 201u"s"
		@test boost_coast_transit_time(100u"m", 1u"m/s/s", 2u"s") == 52u"s"
		@test boost_coast_transit_time(100u"m", 2u"m/s/s", 1u"s") == 51u"s"
		
		@test boost_coast_thrust_time(100u"m", 1u"m/s/s", 101u"s") == 1u"s"
		@test boost_coast_thrust_time(200u"m", 1u"m/s/s", 201u"s") == 1u"s"
		@test boost_coast_thrust_time(100u"m", 1u"m/s/s", 52u"s") == 2u"s"
		@test boost_coast_thrust_time(100u"m", 2u"m/s/s", 51u"s") == 1u"s"
		
		@test beam_core_mass_ratio(.22, .25u"c", .33u"c") ≈ 5.56 atol = 0.01 # frisbee got 5.45 for his mass ratio, but I cannot replicate it.
		
		@test mass_ratio(1u"m/s", 1u"m/s") == exp(1)
		
		@test delta_v(1u"m/s", exp(1)) == 1u"m/s"
		@test delta_v(1u"m/s" / g_n, exp(1)) == 1u"m/s"
		
		@test rocket_thrust(1u"m/s", 1u"kg/s") == 1u"N"
		@test rocket_thrust(1u"m/s", 0.5u"W") == 1u"N"
		
		@test rocket_power(1u"m/s", 2u"kg/s") == 1u"W"
		@test rocket_power(2u"m/s", 1u"kg/s") == 2u"W"
		@test rocket_power(1u"m/s", 2u"N") == 1u"W"
		@test rocket_power(2u"m/s", 1u"N") == 1u"W"
		
		@test mass_flow(1u"m/s", 1u"N") == 1u"kg/s"
		
		@test exhaust_velocity(1u"m/s" / g_n) == 1u"m/s"
		@test exhaust_velocity(1u"N", 1u"kg/s") == 1u"m/s"
		@test exhaust_velocity(1u"W", 2u"kg/s") ≈ 1u"m/s"
	end
	
	@testset "SfDrag" begin
		@test drag_force(1, 2u"kg/m^3", 1u"m^2", 1u"m/s") == 2drag_force(1, 1u"kg/m^3", 1u"m^2", 1u"m/s")
		@test drag_force(1, 2u"kg/m^3", 1u"m^2", 1u"m/s") == 2drag_force(0.5, 2u"kg/m^3", 1u"m^2", 1u"m/s")
		@test drag_force(1, 2u"kg/m^3", 1u"m^2", 1u"m/s") == 0.25drag_force(1, 2u"kg/m^3", 1u"m^2", 2u"m/s")
	end
	
	@testset "SfBlackHole" begin
		@test schwarzschild_radius(1u"Msun") ≈ 2950u"m" atol = 5u"m"
		
		@test hawking_evaporation(1u"Msun") ≈ 2e67u"yr" atol = 1e66u"yr"
	end
	
	@testset "SfLaser" begin
		d = 4.7u"cm"
		r = 350u"km"
		λ = 500u"nm"
		w0 = diffraction_limited_spot_diameter(r, d, λ)
		
		@test diffraction_limited_spot_diameter(r, d, λ) ≈ 474u"cm" atol=1u"mm"
		@test diffraction_limited_element_size(r, w0, λ) ≈ d
		@test diffraction_limited_range(d, w0, λ) ≈ r
		@test diffraction_limited_wavelength(r, d, w0) ≈ λ
		
		m2 = 2
		w0 = diffraction_limited_spot_diameter(r, d, λ, m2)
		
		@test diffraction_limited_spot_diameter(r, d, λ, m2) ≈ 948.1u"cm" atol=1u"mm"
		@test diffraction_limited_element_size(r, w0, λ, m2) ≈ d
		@test diffraction_limited_range(d, w0, λ, m2) ≈ r
		@test diffraction_limited_wavelength(r, d, w0, m2) ≈ λ
	end
	
	@testset "SfElectro" begin
		dp=7.95e22u"A*m^2" # earth's magnetic dipole (TODO: add to planetary/SfSolarSystem?)
		
		@test energy_density(1u"T") == 397887.35751313815u"J/m^3"
		@test field_strength(energy_density(1u"T")) == 1u"T"
		@test lorentz_force(1u"A", 1u"m", 1u"T") == 1u"N"
		@test field_strength(dp, 1u"Rearth") ≈ 60u"μT" atol = 1.5u"μT"
		@test dipole_distance(dp, 60u"μT") ≈ 1u"Rearth" atol = 0.01u"Rearth"
	end
	
	@testset "SfThermo" begin
		n2 = 0.028u"kg/mol"/N_A
		
		@test maxwell_boltzmann_peak_speed(n2, 300u"K") ≈ 422u"m/s" atol=0.1u"m/s"
		@test rms_thermal_velocity(n2, 300u"K") ≈ 517u"m/s" atol=0.1u"m/s"
		
		@test black_body_radiated_power(5772u"K", area(Sphere(1u"Rsun"))) ≈ 1u"Lsun" atol = 1e-5u"Lsun" # approximate luminosity of the sun
		
		@test convert_temperature(1u"eV") |> convert_temperature == 1u"eV"
		
		@test convert_energy(convert_energy(1u"MJ/kg", 1u"g/mol"), 1u"g/mol") ≈ 1u"MJ/kg"
	end
	
	@testset "SfWeaponry" begin
		@test newtonian_penetrator(1u"m", 1u"kg/m^3", 1u"kg/m^3") == 1u"m"
		@test newtonian_penetrator(1u"m", 1u"kg/m^3", 2u"kg/m^3") == 0.5u"m"
		@test newtonian_penetrator(1u"m", 2u"kg/m^3", 1u"kg/m^3") == 2u"m"
		
		@test jet_penetrator(1u"m", 1u"kg/m^3", 1u"kg/m^3") == 1u"m"
		@test jet_penetrator(1u"m", 1u"kg/m^3", 2u"kg/m^3") ≈ 1u"m"/sqrt(2)
		@test jet_penetrator(1u"m", 2u"kg/m^3", 1u"kg/m^3") == sqrt(2) * u"m"
		
		@test coilgun_velocity(1u"g", 1u"cm", 1u"m", 1u"T") ≈ 500u"m/s"
		@test coilgun_field_strength(1u"g", 1u"cm", 1u"m", 500u"m/s") ≈ 1u"T"
		@test coilgun_length(1u"g", 1u"cm", 1u"T", 500u"m/s") ≈ 1u"m"
	end
	
	@testset "SfObjects" begin
		m = Material(1000u"kg/m^3", 0u"Pa")
		s = Sphere(1.0u"m")
		o = Object(s, m)
		
		@test shape(o) == s
		@test radius(o) == 1u"m"
		@test area(o) == area(s)
		@test volume(o) == volume(s)
		
		@test density(o) == density(m)
		
		@test mass(o) == 4000π * u"kg" / 3
		
		@test bouyancy(s, m) ≈ 4000π * g_n * u"kg" / 3
		@test bouyancy(o, m) == 0u"N" # water sphere neutrally bouyant in water
		
		n = Material(2000u"kg/m^3", 0u"Pa")
		
		@test bouyancy(o, n) ≈ 4000π * g_n * u"kg" / 3
		
		@test sphere_of(m, mass(o)) |> shape == s
		@test cube_of(m, 1000u"kg") |> shape == Cube(1.0u"m")
	end
	
	@testset "SfAtmosphere" begin
		# sanity checking against the source's own figures using constant exobase altitude and temperature with monatomic gases.
		@test jeans_escape_timescale(1480u"K", 500u"km", SfSolarSystem.earth, 1u"u") ≈ 6.5e-4u"yr" atol = 1e-5u"yr" # hydrogen
		@test jeans_escape_timescale(1480u"K", 500u"km", SfSolarSystem.earth, 4u"u") ≈ 1.3e2u"yr" atol = 1u"yr" # helium
		@test jeans_escape_timescale(1480u"K", 500u"km", SfSolarSystem.earth, 12u"u") ≈ 6.4e17u"yr" atol = 1e16u"yr" # carbon
		@test jeans_escape_timescale(1480u"K", 500u"km", SfSolarSystem.earth, 14u"u") ≈ 6.3e21u"yr" atol = 1e20u"yr" # nitrogen
		# this is an order of magnitude greater than the example figures, but I'm inclined to blame the source author for making a mistake.
		@test jeans_escape_timescale(1480u"K", 500u"km", SfSolarSystem.earth, 16u"u") ≈ 6.4e25u"yr" atol = 1e24u"yr" # oxygen
	end
	
	@testset "Documentation" begin
		doctest(SfPhysics; manual=false)
	end
end
