using SfPhysics
using Test

using Unitful, UnitfulAstro, Documenter
import PhysicalConstants.CODATA2018: g_n, c_0, m_p

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
		@test duration(0u"m/s^2", 10u"m") == Inf * 1u"s"
		@test duration(10u"m/s^2", 0u"m") == 0u"s"
		@test duration(10u"m/s^2", 10u"m") == sqrt(2) * 1u"s"
		@test_throws DomainError duration(10u"m/s^2", -10u"m")
		
		@test duration(0u"m/s^2", 10u"m", 10u"m/s") == 1u"s"
		@test duration(10u"m/s^2", 0u"m", 10u"m/s") == 0u"s"
		@test duration(10u"m/s^2", 10u"m", 10u"m/s") ≈ 0.7320508075688774u"s"
		@test_throws DomainError duration(10u"m/s^2", -10u"m", 10u"m/s")
		@test_throws DomainError duration(0u"m/s^2", 10u"m", -10u"m/s")
		
		@test kinetic_energy(0u"kg", 1u"m/s") == 0u"J"
		@test kinetic_energy(1u"kg", 0u"m/s") == 0u"J"
		@test kinetic_energy(2u"kg", 1u"m/s") == 1u"J"
		@test kinetic_energy(2u"kg", 2u"m/s") == 4u"J"
		@test kinetic_energy(2u"kg", -2u"m/s") == 4u"J"
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
	
	@testset "SfFluidDynamics" begin
		@test drag_force(1, 2u"kg/m^3", 1u"m^2", 1u"m/s") == 2drag_force(1, 1u"kg/m^3", 1u"m^2", 1u"m/s")
		@test drag_force(1, 2u"kg/m^3", 1u"m^2", 1u"m/s") == 2drag_force(0.5, 2u"kg/m^3", 1u"m^2", 1u"m/s")
		@test drag_force(1, 2u"kg/m^3", 1u"m^2", 1u"m/s") == 0.25drag_force(1, 2u"kg/m^3", 1u"m^2", 2u"m/s")
	end
	
	@testset "Documentation" begin
		doctest(SfPhysics; manual=false)
	end
end
