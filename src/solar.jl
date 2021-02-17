module SfSolarSystem
	
	using Unitful
	using UnitfulAstro
	using UnitfulAngles
	
	using ..SfPlanetary
		
	const sol = Body("Sol", 1u"Msun", 1u"Rsun")

	const mercury = Body(
				"Mercury", 3.3011u"kg" * big(10)^23, 2439.7u"km", 2439.7u"km", 0.088, 
				Orbit(sol, 0.387098u"AU", 0.205630, 174.796, 3.38, 48.331, 23.124),
				Rotation(0.346, 1407.5u"hr", 0.034))

	const venus = Body(
				"Venus", 4.8675u"kg" * big(10)^24, 6051.8u"km", 6051.8u"km", 0.76,
				Orbit(sol, 0.723332u"AU", 0.006772, 50.115, 3.86, 76.680, 54.884),
				Rotation(nothing, -324.025u"d", 177.36))

	const earth = Body(
				"Earth", 1u"Mearth", 1u"Rearth_e", 1u"Rearth_p", 0.306,
				Orbit(sol, 1.000001018u"AU", 0.0167086, 358.617, 7.155, -11.26064, 114.20783), # inclination to ecliptic
				Rotation(0.3307, 0.99726968u"d", 23.4392811))
	const moon = Body(
				"Moon", 7.34767309u"kg" * big(10)^22, 1737.1u"km", 1737.1u"km", 0.136, # not clear if this is bond albedo
				Orbit(earth, 384399u"km", 0.0549,  5.145u"°"),
				Rotation(0.3929, 27.321661u"d", 24))

	const mars = Body(
				"Mars", 6.4171u"kg" * big(10)^23, 3396.2u"km", 3376.2u"km", 0.25,
				Orbit(sol, 227939200u"km", 0.0934, 19.412, 5.65, 49.558, 286.502),
				Rotation(0.3644, 1.025957u"d", 25.19))
	const phobos = Body(
				"Phobos", 1.0659u"kg" * big(10)^16, 11.2667u"km", 11.2667u"km", 0.071, # not clear if this is bond albedo
				Orbit(mars, 9376u"km", 0.0151, 1.093u"°"),
				Rotation(nothing, 0.31891023u"d", 0))
	const deimos = Body(
				"Deimos", 1.4762u"kg" * big(10)^15, 6.2u"km", 6.2u"km", 0.068, # not clear if this is bond albedo
				Orbit(mars, 23463.2u"km", 0.00033, .93u"°"),
				Rotation(nothing, 1.263u"d", 0))

	const ceres = Body(
				"Ceres", 9.3835u"kg" * big(10)^20, 469.73u"km", 469.73u"km", 0.09, # this is geometric albedo; bond albedo not available
				Orbit(sol, 2.7691651545u"AU", 0.0760090291, 77.37209589, 10.59406704, 80.3055316, 73.5976941), # inclination to ecliptic
				Rotation(0.36, 9.07417u"hr", 4))

	const jupiter = Body(
				"Jupiter", 1u"Mjup", 1u"Rjup_e", 1u"Rjup_p", 0.503,
				Orbit(sol, 5.2044u"AU", 0.0489, 20.02, 6.09, 100.464, 273.867),
				Rotation(0.2756, 9.925u"hr", 3.13)) # tilt to orbit

	const saturn = Body(
				"Saturn", 5.6834u"kg" * big(10)^26, 60268u"km", 54364u"km", 0.342,
				Orbit(sol, 9.5826u"AU", 0.0565, 317.02, 5.51, 113.665, 339.392),
				Rotation(0.22, 38018u"s", 26.73)) # tilt to orbit)

	const uranus = Body(
				"Uranus", 8.681u"kg" * big(10)^25, 25559u"km", 24973u"km", 0.3,
				Orbit(sol, 19.2184u"AU", 0.046381, 142.2386, 6.48, 74.006, 96.998857),
				Rotation(0.23, -0.71833u"d", 97.77)) # tilt to orbit

	const neptune = Body(
				"Neptune", 1.02413u"kg" * big(10)^26, 24764u"km", 24341u"km", 0.29,
				Orbit(sol, 30.07u"AU", 0.008678, 256.228, 6.43, 131.784, 276.336),
				Rotation(0.23, 0.6713u"d", 28.32)) # tilt to orbit

	const pluto = Body(
				"Pluto", 1.303u"kg" * big(10)^22, 1188.3u"km", 1188.3u"km", 0.575, # average geometric; bond albedo not available
				Orbit(sol, 39.482u"AU", 0.2488, 14.53, 11.88, 110.299, 113.834),
				Rotation(nothing, 6.3872304u"d", 122.53)) # tilt to orbit
	const charon = Body(
				"Charon", 1.586u"kg" * big(10)^21, 606u"km", 606u"km", 0.35, # average geometric; bond albedo not available
				Orbit(pluto, 19591.4u"km", 0.0002, nothing, 0.08, 223.046, nothing), # SMA to centre of Pluto, not to common barycenter
				Rotation(nothing, 6.3872304u"d", 0))
			
end