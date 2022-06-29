# Conversion between Seis types

# Can only convert to the same subtype of Position
for Geom in (Cartesian, Geographic)
    geom = nameof(Geom)
    @eval begin
        # Identity 'conversions'
        Base.convert(::Type{$geom{T}}, p::$geom{T}) where T = p
        Base.convert(::Type{Event{T,$geom{T}}}, e::Event{T,$geom{T}}) where T = e
        Base.convert(::Type{Station{T,$geom{T}}}, s::Station{T,$geom{T}}) where T = s

        # Parameter T conversions
        Base.convert(::Type{$geom{Tout}}, p::$geom{Tin}) where {Tin,Tout} =
            $geom{Tout}($([:(getfield(p, $i)) for i in 1:fieldcount(Geom)]...))

        Base.convert(::Type{Event{Tout,$geom{Tout}}}, e::Event{Tin,$geom{Tin}}) where {Tin, Tout} =
            Event{Tout,$geom{Tout}}(e.pos, e.time, e.id, e.meta)

        Base.convert(::Type{Station{Tout,$geom{Tout}}}, s::Station{Tin,$geom{Tin}}) where {Tin, Tout} =
            Station{Tout,$geom{Tout}}(s.net, s.sta, s.loc, s.cha, s.pos, s.elev,
                s.azi, s.inc, s.meta)

        Base.convert(::Type{Trace{T,V,$geom{T}}}, t::Trace{T,V,$geom{T}}) where {T,V} = t
        function Base.convert(::Type{Trace{Tout,Vout,$geom{Tout}}},
                t::Trace{Tin,Vin,$geom{Tin}}
                ) where {Tin,Vin,Tout,Vout}
            Trace{Tout,Vout,$geom{Tout}}($([:(getfield(t, $i)) for i in 1:fieldcount(Trace)]...))
        end

        Base.convert(::Type{FourierTrace{T,V,$geom{T}}}, f::FourierTrace{T,V,$geom{T}}) where {T,V} = f
        function Base.convert(::Type{FourierTrace{Tout,Vout,$geom{Tout}}},
                f::FourierTrace{Tin,Vin,$geom{Tin}}
                ) where {Tin,Tout,Vin,Vout}
            FourierTrace{Tout,Vout,$geom{Tout}}($([:(getfield(f, $i)) for i in 1:fieldcount(FourierTrace)]...))
        end
    end
end
