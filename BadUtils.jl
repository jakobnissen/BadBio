module BadUtils

using DataStructures

nested_dict(K::DataType=Any, V::DataType=Any) = DefaultDict{K, V}(nested_dict)

function reversekv(dict::AbstractDict{K,V}) where {K,V}
    Base.typename(typeof(dict)).wrapper{V,K}(v => k for (k, v) in dict)
end

collectkeys(args) = keys(args) |> collect
collectvals(args) = values(args) |> collect

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

function tryint(number)
    return (try
        Int64(number)
    catch
        number
    end)
end
function tryfloat(number)
    return (try
        Float64(number)
    catch
        number
    end)
end
function tryfloat32(number)
    return (try
        Float32(number)
    catch
        number
    end)
end
function unitrange(arr::AbstractVector{T}) where {T<:Int}
    start_value = arr[1]
    end_value = arr[end]
    end_value == start_value && error("the start and end points are the same value")
	end_value < start_value && error("the last value is lower than the first")
    for i in 2:size(arr,1)
        if arr[i]-arr[i-1] != 1
            error("inconsistent step for unit range")
        end
    end
    return UnitRange(start_value,end_value)
end
function steprange(arr::AbstractVector{T}) where {T<:Real}
    start_value = arr[1]
    end_value = arr[end]
	start_value == end_value && error("the start and end points are the same value")
    step = 1
    for i in 2:size(arr,1)
        if i == 2
            step = arr[i]-arr[i-1]
        end
        if arr[i]-arr[i-1] != step
            error("inconsistent step for step range")
        end
    end
    return StepRangeLen(start_value,step,length(arr))
end
function varcall(name::String, body::Any)
    name=Symbol(name)
    @eval (($name) = ($body))
	return Symbol(name)
end

end # module
