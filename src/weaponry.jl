using Unitful

function cratering_strength(yield_strength::Unitful.Pressure)
    3u"J/m^3" * Unitful.ustrip(u"Pa", yield_strength)
end
function cratering_volume(e::Unitful.Energy, yield_strength::Unitful.Pressure)
    e / cratering_strength(yield_strength)
end
function cratering_depth(e::Unitful.Energy, yield_strength::Unitful.Pressure)
    r = cratering_volume(e, yield_strength)
    # hemisphere volume = (2/3)πr^3

    cbrt(3r/2π)
end

jet_penetrator(l::Unitful.Length, ρj::Unitful.Density, ρt::Unitful.Density) = l * sqrt(ρj / ρt)
