module SfObjects

using Unitful, ..SfGeometry, ..SfMatter

import ..SfMatter: mass, density

export mass

mass(m::Material, s::Shape) = density(m) * volume(s) |> u"kg"
mass(d::Unitful.Density, s::Shape) = d * volume(s) |> u"kg"

end
