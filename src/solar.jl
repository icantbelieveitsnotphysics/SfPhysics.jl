module SfSolarSystem
	
using Unitful, UnitfulAstro, UnitfulAngles, Measurements

using ..SfPlanetary, ..SfGeometry
	
# useful sources:
# - mostly wikipedia
# - https://ssd.jpl.nasa.gov/?sat_phys_par (planetary satellite physical parameters)
# - https://arxiv.org/abs/1604.06129 (pluto and charon albedo)
	
const sol = Star("Sol", 1u"Msun", 1u"Rsun", 5772u"K", 4.83, "G2V")

const mercury = Body(
			"Mercury", (0.330114±0.000021)*1e24u"kg", Sphere((2439.4±0.1)u"km"), 0.088, 0.142,
			Orbit(sol, 0.387098u"AU", 0.205630, 3.38, 48.331, 23.124),
			Rotation(0.346, 1407.5u"hr", 0.034))

const venus = Body(
			"Venus", 4.8675u"kg" * big(10)^24, Sphere(6051.8u"km"), 0.76, 0.689,
			Orbit(sol, 0.723332u"AU", 0.006772, 3.86, 76.680, 54.884),
			Rotation(missing, -324.025u"d", 177.36))

const earth = Body(
			"Earth", 1u"Mearth", Spheroid(1u"Rearth_e", 1u"Rearth_p"), 0.306, 0.367,
			Orbit(sol, 1.000001018u"AU", 0.0167086, 7.155, -11.26064, 114.20783), # inclination to ecliptic
			Rotation(0.3307, 0.99726968u"d", 23.4392811))
const moon = Body(
			"Moon", 7.34767309u"kg" * big(10)^22, Sphere(1737.1u"km"), 0.11, 0.12,
			Orbit(earth, 384399u"km", 0.0549,  5.145u"°"),
			Rotation(0.3929, 27.321661u"d", 24))

const mars = Body(
			"Mars", 6.4171u"kg" * big(10)^23, Spheroid(3396.2u"km", 3376.2u"km"), 0.25, 0.17,
			Orbit(sol, 227939200u"km", 0.0934, 5.65, 49.558, 286.502),
			Rotation(0.3644, 1.025957u"d", 25.19))
const phobos = Body(
			"Phobos", 1.0659u"kg" * big(10)^16, Ellipsoid(13.5u"km", 11u"km", 9u"km"), missing, 0.071±0.012,
			Orbit(mars, 9376u"km", 0.0151, 1.093u"°"),
			Rotation(missing, 0.31891023u"d", 0))
const deimos = Body(
			"Deimos", 1.4762u"kg" * big(10)^15, Ellipsoid(7.5u"km", 6.1u"km", 5.5u"km"), missing, 0.068±0.007,
			Orbit(mars, 23463.2u"km", 0.00033, .93u"°"),
			Rotation(missing, 1.263u"d", 0))

const ceres = Body(
			"Ceres", 9.3835u"kg" * big(10)^20, Sphere(469.73u"km"), missing, 0.090±0.0033,
			Orbit(sol, 2.7691651545u"AU", 0.0760090291, 10.59406704, 80.3055316, 73.5976941), # inclination to ecliptic
			Rotation(0.36, 9.07417u"hr", 4))

const jupiter = Body(
			"Jupiter", 1u"Mjup", Spheroid(1u"Rjup_e", 1u"Rjup_p"), 0.503, 0.538,
			Orbit(sol, 5.2044u"AU", 0.0489, 6.09, 100.464, 273.867),
			Rotation(0.2756, 9.925u"hr", 3.13)) # tilt to orbit
			
const io = Body(
			"Io", 8.931938e22u"kg", Sphere(18621.6u"km"), missing, 0.63±0.02,
			Orbit(jupiter, 421700u"km", 0.0041, 0.05), # not quite an elliptical orbit; has resonance issues
			Rotation(0.37824, 1.769137786u"d", 0))

const europa = Body(
			"Europa", 4.799844e22u"kg", Sphere(1560.8u"km"), missing, 0.67±0.03,
			Orbit(jupiter, 670900u"km", 0.009, 0.47), # not quite an elliptical orbit
			Rotation(0.346, 3.551181u"d", 0.1))

const ganymede = Body(
			"Ganymede", 1.4819e23u"kg", Sphere(2634.1u"km"), missing, 0.43±0.02,
			Orbit(jupiter, 1070400u"km", 0.0013, 0.2), # this one is elliptical
			Rotation(0.3115, 7.15455296u"d", 0.33)),

const callisto = Body(
			"Callisto", 1.075938e23u"kg", Sphere(2410.3u"km"), missing, 0.17±0.02,
			Orbit(jupiter, 1882700u"km", 0.0074, 0.192), # also actually elliptical
			Rotation(0.3549, 16.6890184u"d", 0))

const saturn = Body(
			"Saturn", 5.6834u"kg" * big(10)^26, Spheroid(60268u"km", 54364u"km"), 0.342, 0.499,
			Orbit(sol, 9.5826u"AU", 0.0565, 5.51, 113.665, 339.392),
			Rotation(0.22, 38018u"s", 26.73)) # tilt to orbit)
			
const titan = Body(
			"Titan", 1.3452e23u"kg", Sphere(2574.73u"km"), missing, 0.22,
			Orbit(saturn, 1221870u"km", 0.0288, 0.34854),# inclination to Saturn's equator
			Rotation(0.3414, 15.945u"d", 0))

const uranus = Body(
			"Uranus", 8.681u"kg" * big(10)^25, Spheroid(25559u"km", 24973u"km"), 0.3, 0.488,
			Orbit(sol, 19.2184u"AU", 0.046381, 6.48, 74.006, 96.998857),
			Rotation(0.23, -0.71833u"d", 97.77)) # tilt to orbit

const neptune = Body(
			"Neptune", 1.02413u"kg" * big(10)^26, Spheroid(24764u"km", 24341u"km"), 0.29, 0.422,
			Orbit(sol, 30.07u"AU", 0.008678, 6.43, 131.784, 276.336),
			Rotation(0.23, 0.6713u"d", 28.32)) # tilt to orbit

const pluto = Body(
			"Pluto", 1.303u"kg" * big(10)^22, Sphere(1188.3u"km"), 0.72±0.07, 0.575, # average geometric
			Orbit(sol, 39.482u"AU", 0.2488, 11.88, 110.299, 113.834),
			Rotation(missing, 6.3872304u"d", 122.53)) # tilt to orbit
const charon = Body(
			"Charon", 1.586u"kg" * big(10)^21, Sphere(606u"km"), 0.25±0.03 , 0.35, # average geometric
			Orbit(pluto, 19591.4u"km", 0.0002, 0.08, 223.046, missing), # SMA to centre of Pluto, not to common barycenter
			Rotation(missing, 6.3872304u"d", 0))
			
end
