"""
    SeisDict

Wrapper around `Base.Dict` which allows one to get or set values using
the `{get|set}property[!]` syntax, i.e., `.`-access like `dict.key`.
`SeisDict`s also differ in that `missing` is returned instead of throwing a
`KeyError` when accessing a nonexistent key.  A key is removed if it is set to
`missing`.

# Examples

```
julia> dict = Seis.SeisDict(:a=>1)
Seis.SeisDict{Any,Any} with 1 entry:
  :a => 1

julia> dict.a
1

julia> dict.b
missing

julia> dict.a = missing
missing

julia> dict
Seis.SeisDict{Any,Any} with 0 entries
```

Access via `getindex` and `setindex!` (using `[]`s) is still possible:

```
julia> dict[:c] = 3
3

julia> dict[:c]
3

julia> dict[:c] = missing
missing

julia> dict[:d] # No key with this value
missing

julia> dict
Seis.SeisDict{Any,Any} with 0 entries
```

Arrays of `SeisDict`s also define `.`-access and setting via broadcasting, so one
may do:

```
julia> d1 = Seis.SeisDict(:a=>1)
Seis.SeisDict{Any,Any} with 1 entry:
  :a => 1

julia> d2 = deepcopy(d1); d2.a = 2;

julia> d = [d1, d2]
2-element Array{Seis.SeisDict{Any,Any},1}:
 Seis.SeisDict(:a => 1)
 Seis.SeisDict(:a => 2)

julia> d.a
2-element Array{Int64,1}:
 1
 2
```
"""
struct SeisDict{K,V} <: Base.AbstractDict{K,V}
    dict::Dict{K,V}
end

SeisDict{K,V}(args...) where {K,V} = SeisDict{K,V}(Dict{K,V}(args...))
SeisDict{K}(args...) where {K} = SeisDict{K,Any}(Dict{K,Any}(args...))
SeisDict(args...) = SeisDict{Any,Any}(Dict{Any,Any}(args...))

# Have to use getfield here because we define getproperty as getindex
Base.getindex(sd::SeisDict{K,V}, key) where {K,V} = get(getfield(sd, :dict), key, missing)
Base.setindex!(sd::SeisDict{K,V}, ::Missing, key) where {K,V} = delete!(getfield(sd, :dict), key)
Base.setindex!(sd::SeisDict{K,V}, v, i...) where {K,V} = setindex!(getfield(sd, :dict), v, i...)
Base.Dict(sd::SeisDict) = getfield(sd, :dict)
Base.empty!(sd::SeisDict) = empty!(getfield(sd, :dict))
Base.summary(sd::SeisDict) = summary(getfield(sd, :dict))
Base.iterate(sd::SeisDict, args...) = iterate(getfield(sd, :dict), args...)
Base.length(sd::SeisDict) = length(getfield(sd, :dict))
Base.get(sd::SeisDict, args...) = get(getfield(sd, :dict), args...)
Base.delete!(sd::SeisDict, args...) = (delete!(getfield(sd, :dict), args...); sd)

Base.getproperty(sd::SeisDict, key::Symbol) = getindex(sd, key)
Base.setproperty!(sd::SeisDict{K,V}, key::Symbol, ::Missing) where {K,V} = delete!(sd, key)
Base.setproperty!(sd::SeisDict{K,Any}, key::Symbol, ::Missing) where K = delete!(sd, key)
Base.setproperty!(sd::SeisDict{K,V}, key::Symbol, val) where {K,V} = setindex!(sd, V(val), key)
Base.setproperty!(sd::SeisDict{K,Any}, key::Symbol, val) where K = setindex!(sd, val, key)
Base.propertynames(sd::SeisDict, private::Bool=false) = collect(keys(sd))

Base.getproperty(sd::AbstractArray{<:SeisDict}, key::Symbol) = getproperty.(sd, key)
Base.setproperty!(sd::AbstractArray{<:SeisDict}, key::Symbol, val) = setproperty!.(sd, key, val)
