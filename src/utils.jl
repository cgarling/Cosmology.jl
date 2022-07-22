
isscalar(x) = isa(x, Union{Number,AbstractString,Char,Bool})::Bool
""" 
    heaviside(x::Union{AbstractFloat,AbstractArray}; y=0.5)
Returns the heaviside function which is 0 if x<0, 1 if x>0, and y if x==0. """
heaviside(x::AbstractFloat, x2=0.5) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,x2)))
heaviside(x::AbstractArray, x2=0.5) = [heaviside(i,x2) for i in x]

##############################################################################
"""
    logspace(start,stop; length=nothing, step=nothing)
Returns an array for numbers equally spaced in log10 between start and stop. Step is the difference between adjacent values in log10, or length will determine the step on its own.
"""
logspace(start,stop; length=nothing, step=nothing) = exp10.( range(log10(start), log10(stop), step=step, length=length) )


"""
    ln_derivative(f_x,df_dx)
Given the value of a function `f_x` and the point derivative `df_dx` calculate the ln derivative as `df_dx / f_x`. """
ln_derivative(f_x,df_dx) = df_dx ./ f_x

"""
    lnln_derivative(x,f_x,df_dx)
Given the value of a function `f_x` and the point derivative `df_dx` at `x` calculate the dimensionless loglog derivative as `df_dx / (x * f_x)`. """
lnln_derivative(x,f_x,df_dx) = df_dx ./ (x .* f_x)

"""
    cumtrapz(X::AbstractVector, Y::AbstractVector;normalized=false)
Simple function to compute the cumulative integral over a vector of function values `Y` evaluated at positions `X`. If `normalized` is `true`, then returns the cumulative integral normalized so that the final value is 1. Works with Measurements.jl. """
function cumtrapz(X::AbstractVector, Y::AbstractVector;normalized=false) # where {T <: AbstractVector}
    @assert length(X) == length(Y)
    out = similar(X)
    out[1] = zero(eltype(X))
    # Iterate over arrays
    @inbounds for i in 2:length(X)
        out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    normalized ? out./out[end] : out
end

""" 
    smooth_array_gaussian(x::AbstractVector,y::AbstractVector,σ::Real;nsig::Real=3)
Smooth a vector `y` using a Gaussian kernel with standard deviation `σ`. `σ` is defined in the units of `x`. Will not smooth points within `nsig * sigma` units of the ends. """
function smooth_array_gaussian(x::AbstractVector,y::AbstractVector,σ::Real;nsig::Real=3)
    @assert length(x) == length(y)
    result=zero(y)
    @inbounds for i in eachindex(x,y)
        total=zero(eltype(y))
        if abs(x[i]-x[1]) < nsig*σ || abs(x[i]-x[end]) < nsig*σ
            result[i]=y[i]
        else
            @inbounds for j in eachindex(x,y)
                dist = x[i]-x[j]
                abs(dist)>5*σ ? weight=0.0 : weight = exp(-(dist)^2/(2*σ^2))
                # weight = exp(-(x[i]-x[j])^2/(2*σ^2))
                result[i] += y[j] * weight
                total += weight
            end
            result[i] /= total
        end
    end
    return result
end
