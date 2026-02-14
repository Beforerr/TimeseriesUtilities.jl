# Common API for different types of arrays
# AxisKeys.jl: https://github.com/mcabbott/AxisKeys.jl
# DimensionalData.jl: https://github.com/rafaqz/DimensionalData.jl

"""
# TimeseriesUtilities.jl Array Interface

This module defines a generic interface for working with multi-dimensional timeseries arrays.
To support a custom array type, implement the following functions:

## Required Interface Functions

### Dimension Querying
- `dimnum(x, query)`: Get the numeric index of a dimension given a query (e.g., :time or TimeDim)
- `dims(x, dim)`: Get dimension object(s) for the given dimension index/query
- `axiskeys(x, dim)`: Get the coordinate values (axis keys) for a dimension

### Data Access
- `unwrap(x)`: Extract the raw underlying array data (without metadata)
- `unwrap(dim::Dimension)`: Extract coordinate values from a dimension object

### Array Reconstruction
- `rebuild(x, data)`: Rebuild array with new data, preserving dimensions
- `rebuild(x, data, dims)`: Rebuild array with new data and new dimensions
- `rebuild(dim::Dimension, values)`: Rebuild a dimension with new coordinate values
- `set(x, dim => values)`: Update a dimension's coordinate values

### Time-Specific (Optional, defaults provided)
- `times(x)`: Get the time coordinate values (defaults to axiskeys(x, dimnum(x, nothing)))
- `timedim(x, query)`: Get the time dimension object

## Example Extension

```julia
module MyArrayExt
import TimeseriesUtilities: unwrap, dimnum, set, dims, axiskeys, rebuild

unwrap(x::MyArray) = x.data
dimnum(x::MyArray, query) = findfirst(==(query), x.dimnames)
axiskeys(x::MyArray, dim) = x.axes[dim]
dims(x::MyArray, dim) = x.dimnames[dim]

rebuild(x::MyArray, data) = MyArray(data, x.axes, x.dimnames)
rebuild(x::MyArray, data, newdims) = MyArray(data, newdims, x.dimnames)

function set(x::MyArray, pair::Pair)
    dim, newkeys = pair
    dn = dimnum(x, dim)
    newaxes = ntuple(i -> i == dn ? newkeys : x.axes[i], ndims(x))
    return MyArray(x.data, newaxes, x.dimnames)
end
end
```
"""
module API end

"""
    dimnum(x, query)

Get the numeric index of a dimension given a query.

# Arguments
- `x`: Array with dimensions
- `query`: Dimension query (e.g., :time, TimeDim, or nothing for default time dimension)

# Returns
- Integer index of the dimension

# Example
```julia
dimnum(my_array, :time)  # Returns 1 if time is the first dimension
```

Extend this function for custom array types. See module docstring for examples.
"""
function dimnum end

"""
    set(x, dim => values)

Update a dimension's coordinate values in an array.

# Arguments
- `x`: Array with dimensions
- `dim => values`: Pair of dimension query and new coordinate values

# Returns
- New array with updated dimension coordinates

# Example
```julia
set(my_array, :time => new_times)
```

Extend this function for custom array types. See module docstring for examples.
"""
function set end

"""
    axiskeys(x, dim)

Get the coordinate values (axis keys) for a dimension.

# Arguments
- `x`: Array with dimensions
- `dim`: Dimension index or query

# Returns
- Vector of coordinate values for the dimension

# Example
```julia
axiskeys(my_array, 1)  # Get coordinates for first dimension
axiskeys(my_array, :time)  # Get time coordinates
```

Extend this function for custom array types. See module docstring for examples.
"""
function axiskeys end

"""
    dims(x, dim)

Get dimension object(s) for the given dimension index/query.

# Arguments
- `x`: Array with dimensions
- `dim`: Dimension index or query

# Returns
- Dimension object (implementation-specific)

# Example
```julia
dims(my_array, 1)  # Get first dimension
dims(my_array, :time)  # Get time dimension
```

Extend this function for custom array types. See module docstring for examples.
"""
function dims end

"""
    unwrap(x)

Extract the raw underlying array data without metadata.

# Arguments
- `x`: Array with dimensions or dimension object

# Returns
- Raw array data (for arrays) or coordinate values (for dimensions)

# Example
```julia
unwrap(my_array)  # Returns plain Array
unwrap(time_dim)  # Returns time coordinate vector
```

Extend this function for custom array types. See module docstring for examples.
"""
function unwrap end

"""
    rebuild(x, data)
    rebuild(x, data, dims)
    rebuild(dim::Dimension, values)

Rebuild an array or dimension with new data/values.

# Arguments
- `x`: Array with dimensions or dimension object
- `data`: New data array
- `dims`: (Optional) New dimension objects
- `values`: New coordinate values (for dimension rebuild)

# Returns
- New array/dimension with updated data/values, preserving metadata

# Example
```julia
rebuild(my_array, new_data)  # New array with same dimensions
rebuild(my_array, new_data, new_dims)  # New array with new dimensions
rebuild(time_dim, new_times)  # New dimension with new coordinates
```

Extend this function for custom array types. See module docstring for examples.
"""
function rebuild end

"""
    times(x, args...)

Get the time coordinate values from an array.

# Arguments
- `x`: Array with time dimension
- `args...`: Additional arguments passed to `timedim`

# Returns
- Vector of time coordinates

# Example
```julia
times(my_array)  # Get default time coordinates
times(my_array, :timestamp)  # Get specific time dimension
```

This function has a default implementation but can be overridden for custom types.
"""
function times end

"""
    timedim(x, query=nothing)

Get the time dimension object from an array.

# Arguments
- `x`: Array with time dimension
- `query`: (Optional) Dimension query for time (defaults to implementation-specific time dimension)

# Returns
- Time dimension object

# Example
```julia
timedim(my_array)  # Get default time dimension
timedim(my_array, :timestamp)  # Get specific time dimension
```

Extend this function for custom array types. See module docstring for examples.
"""
function timedim end
