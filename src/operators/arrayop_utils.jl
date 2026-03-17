"""
Utility functions for constructing and evaluating ArrayOp expressions.
"""

const ConstSR = BSImpl.Const{SymReal}

"""
    const_wrap(arr)

Wrap a numerical array in `Const{SymReal}` for use in ArrayOp indexing.
"""
const_wrap(arr) = ConstSR(arr)

"""
    get_idx_vars(ndim)

Get `ndim` symbolic index variables from the shared ArrayOp index pool.
"""
function get_idx_vars(ndim::Int)
    idxs_arr = idxs_for_arrayop(SymReal)
    return [idxs_arr[d] for d in 1:ndim]
end

"""
    make_arrayop(idx_vars, expr, ranges)

Construct an `ArrayOp{SymReal}`.
"""
function make_arrayop(idx_vars, expr, ranges)
    return SymbolicUtils.ArrayOp{SymReal}(idx_vars, unwrap(expr), +, nothing, ranges)
end

"""
    evaluate_arrayop(ao)

Evaluate an ArrayOp built from `Const`-wrapped data to a numerical array.

Scalarizes the ArrayOp and extracts Float64 values. For testing and validation.
"""
function evaluate_arrayop(ao)
    s = Symbolics.scalarize(wrap(ao))
    return Float64.(Symbolics.value.(s))
end
