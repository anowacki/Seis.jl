# Define hash so that Base.isequal (and so Base.unique) work correctly
# for our types.  SeisDict does the right thing already.
for T in (Cartesian, Geographic, Event, Station, Trace)
    Tsym = nameof(T)
    fields = fieldnames(T)
    @eval function Base.hash(x::$Tsym, h::UInt)
        $([:(h = hash(x.$f, h)) for f in fields]...)
        h
    end
end
