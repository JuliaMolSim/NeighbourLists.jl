# Adapt.adapt_structure rules for the public neighbour-list types.
#
# Without these rules, passing a `SortedCellList` or `PairList` as an
# argument to a `KernelAbstractions.@kernel` fails GPU compilation with
# `KernelError: passing non-bitstype argument` — CUDA.jl's
# `KernelAdaptor` doesn't recursively adapt the inner CuArray fields to
# their `CuDeviceArray` counterparts without an explicit
# `adapt_structure` method on the wrapping struct.
#
# Internal NL.jl kernels (in `gpu_kernels.jl`) avoid this by unpacking
# the struct into separate kernel arguments. External consumers calling
# `for_each_neighbour(clist, i)` inside a custom `@kernel` need the
# struct itself to be device-friendly. The rules below recursively
# adapt the array fields while keeping the immutable scalar / SVector
# fields as-is.
#
# See issue #37.

using Adapt

Adapt.adapt_structure(to, clist::SortedCellList) =
    SortedCellList(
        Adapt.adapt(to, clist.X),
        Adapt.adapt(to, clist.X_orig),
        Adapt.adapt(to, clist.perm),
        Adapt.adapt(to, clist.cell_id),
        Adapt.adapt(to, clist.cell_offsets),
        clist.cell,
        clist.inv_cell,
        clist.pbc,
        clist.cutoff,
        clist.ncells,
        clist.ncells_total,
    )

Adapt.adapt_structure(to, p::PairList) =
    PairList(
        Adapt.adapt(to, p.X),
        p.C,
        p.cutoff,
        Adapt.adapt(to, p.i),
        Adapt.adapt(to, p.j),
        Adapt.adapt(to, p.S),
        Adapt.adapt(to, p.first),
    )
