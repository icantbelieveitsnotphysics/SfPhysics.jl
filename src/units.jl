module SfUnits

using Unitful

export kardashev, tnt

@dimension KP "KP" KardashevPower
@derived_dimension KPower KP

const Power{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐋^2*Unitful.𝐌*Unitful.𝐓^-3,U}
const Energy{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐋^2*Unitful.𝐌*Unitful.𝐓^-2,U}
const SpecificEnergy{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐋^2*Unitful.𝐓^-2,U}
const AngularVelocity{T, U} = Unitful.AbstractQuantity{T, Unitful.𝐓 ^-1, U}
const AngularMomentum{T, U} = Unitful.AbstractQuantity{T, Unitful.𝐋^2 * Unitful.𝐌 * Unitful.𝐓 ^-1, U}
const Acceleration{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐋*Unitful.𝐓^-2, U}
const Speed{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐋*Unitful.𝐓^-1, U}
const Mass{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐌, U}
const MolarMass{T, U} = Unitful.AbstractQuantity{T,Unitful.𝐌 * Unitful.𝐍^-1, U}

const Angle{T} = Union{ Quantity{T, NoDims, typeof(u"°")}, Quantity{T, NoDims, typeof(u"rad")} }
	
to_angle(a::Real) = Angle(a * u"°")	
to_angle(a::Angle) = a
to_angle(::Nothing) = nothing
to_angle(::Missing) = missing

@unit tt "tt" TonneTNT 4.184e+9u"J" true
@refunit Kdv "Kdv" Kardashev KP false

function __init__()
    Unitful.register(SfUnits)
end

Unitful.register(SfUnits)

"""
    kardashev(p::Power)
	
Convert a power `p` (eg. in watts) to its equivalent in the Kardashev scale.
"""
kardashev(p::Power) =  u"Kdv" * (log10(ustrip(p |> u"W")) - 7.38) / 9.73

"""
    kardashev(k::KPower)

Convert a Kardashev scale value into watts.
"""
kardashev(k::KPower) = u"W" * 10^((ustrip(k |> u"Kdv") * 9.73) + 7.38)

"""
    tnt(e::Energy)
	
Convenience method for converting an energy quantity into tonnes of TNT equivalent.
"""
tnt(e::Energy) = e |> u"tt"

end
