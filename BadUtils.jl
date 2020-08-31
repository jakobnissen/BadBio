module BadUtils

using DataStructures

nested_dict(K::DataType=Any, V::DataType=Any) = DefaultDict{K, V}(nested_dict)

"""
    Path(s::String)

Create a path from a `String`. Paths may be concatenated with `*`
"""
struct Path
    x::String
end

Base.convert(::Type{String}, p::Path) = p.x
Base.String(p::Path) = convert(String, p)
Base.open(x::Path, m::String) = open(convert(String, x), m)

# Necessary to specify Function due to ambiguity errors
Base.open(f::Function, x::Path) = open(f, convert(String, x))
Base.:*(x::Path, y::Union{Path, String}) = Path(joinpath(String(x), String(y)))
Base.split(x::Path) = splitpath(x.x)

for f in [:(Base.open), :(Base.splitpath)]
    @eval $f(x::Path) = $f(String(x))
end

for f in [:(Base.basename), :(Base.isdir), :(Base.ispath), :(Base.isfile)]
    @eval $f(x::Path) = Path($f(String(x)))
end

"""
    @indir directory expression

Move to `directory`, execute expression, then move back.
"""
macro indir(directory, ex)
    quote
        local dir = pwd()
        cd($(directory))
        local val = $(esc(ex))
        cd(dir)
        val
    end
end

end # module
