
using Base.Threads

export maptosites!, maptosites_d!
export map_sites!, map_pairs!, map_pairs_d!, map_pairs_d_vec!


function mt_split(niter::TI, maxthreads=MAX_THREADS[1]) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   # nn = ceil.(TI, linspace(1, niter+1, nt+1))
   nn = ceil.(TI, range(1, stop=niter+1, length=nt+1))
   rgs = [nn[i]:(nn[i+1]-1) for i = 1:nt]
   return nt, rgs
end

function mt_split_interlaced(niter::TI, maxthreads=MAX_THREADS[1]) where TI
   nt = minimum([maxthreads, nthreads(), niter])
   rgs = [ j:nt:niter for j = 1:nt ]
   return nt, rgs
end


function _mt_map_!(f::FT, out, it, inner_loop) where FT
   nt, rg = mt_split(length(it))
   if nt == 1
      inner_loop(f, out, it, 1:length(it))
   else
      OUT = [[out]; [zeros(out) for i = 2:nt]]
      @threads for i = 1:nt
         inner_loop(f, OUT[i], it, rg[i])
      end
      for it = 2:nt
         out .+= OUT[it]
      end
   end
   return out
end

maptosites!(f, out::AbstractVector, it::AbstractIterator) =
   _mt_map_!(f, out, it, maptosites_inner!)

maptosites_d!(f, out::AbstractVector, it::AbstractIterator) =
   _mt_map_!(f, out, it, maptosites_d_inner!)

# ============ assembly over pairs


"""
`mapreduce_sym!{S, T}(out::AbstractVector{S}, f, it::PairIterator{T})`

symmetric variant of `mapreduce!{S, T}(out::AbstractVector{S}, ...)`, summing only
over bonds (i,j) with i < j and adding f(R_ij) to both sites i, j.
"""
function maptosites_inner!(f::FT, out, it::PairIterator, rg) where FT
   nlist = it.nlist
   for n in rg
      if  nlist.i[n] < nlist.j[n]
         f_ = f(nlist.r[n], nlist.R[n]) / 2
         out[nlist.i[n]] += f_
         out[nlist.j[n]] += f_
      end
   end
   return out
end


"""
`mapreduce_antisym!{T}(out::AbstractVector{SVec{T}}, df, it::PairIterator{T})`

anti-symmetric variant of `mapreduce!{S, T}(out::AbstractVector{S}, ...)`, summing only
over bonds (i,j) with i < j and adding f(R_ij) to site j and
-f(R_ij) to site i.
"""
function maptosites_d_inner!(f::FT, out, it::PairIterator, rg) where FT
   nlist = it.nlist
   for n in rg
      if nlist.i[n] < nlist.j[n]
         f_ = f(nlist.r[n], nlist.R[n])
         out[nlist.j[n]] += f_
         out[nlist.i[n]] -= f_
      end
   end
   return out
end


# ============ assembly over sites

function maptosites!(f::FT, out::AbstractVector, it::SiteIterator) where FT
   @threads for i = 1:nsites(it.nlist)
      _, R = neigs(it.nlist, i)
      out[i] = f(R)
   end
   return out
end

function maptosites_d!(df::FT, out::AbstractVector, it::SiteIterator) where FT
   nt = nthreads()
   OUT = [out; [zeros(out) for n = 2:nt]]
   @threads for i = 1:nsites(it.nlist)
      j, R = neigs(it.nlist, i)
      df_ = df(R)
      OUT[threadid()][j] += df_
      OUT[threadid()][i] -= sum(df_)
   end
   for it = 2:nt  # Fixed: was `n`, should be `nt`
      out .+= OUT[it]
   end
   return out
end


# ==================== Threading Helper ====================

"""
    threaded_accumulate!(f, out, range)

Execute `f(thread_out, item)` for each item in `range` in parallel,
accumulating results into `out` using thread-local copies.

`f(thread_out, item)` should modify `thread_out` in-place for each item.
Results from all threads are reduced back into `out`.
"""
function threaded_accumulate!(f::F, out::AbstractVector, range) where F
    nt = nthreads()
    if nt == 1
        for item in range
            f(out, item)
        end
    else
        max_tid = Threads.maxthreadid()
        outs = [i == 1 ? out : zeros(eltype(out), length(out)) for i in 1:max_tid]
        @threads for item in range
            tid = threadid()
            f(outs[tid], item)
        end
        for tid in 2:max_tid
            out .+= outs[tid]
        end
    end
    return out
end


# ==================== KernelAbstractions-based operations ====================


"""
    map_sites!(f, out, clist::SortedCellList)

Apply function `f` to each site's neighbourhood and store in `out`.
`f(Rs)` receives a vector of displacement vectors to neighbours.

Uses KernelAbstractions for parallelization.
"""
function map_sites!(f::F, out::AbstractVector, clist::SortedCellList{T,TI}) where {F, T, TI}
    backend = get_array_backend(out)
    nat = nsites(clist)

    # For CPU backend, use threaded loop
    if backend isa CPU
        @threads for i in 1:nat
            js, Rs, Ss = neighbours(clist, i)
            out[i] = f(Rs)
        end
    else
        # GPU kernel not feasible: f(Rs) expects dynamically-sized vector,
        # which is incompatible with GPU execution. For GPU workflows,
        # use map_pairs! or map_pairs_d! instead.
        error("map_sites! is not supported on GPU. Use map_pairs! or map_pairs_d! for GPU-accelerated operations.")
    end

    return out
end

"""
    map_sites!(f, out, nlist::PairList)

Apply function `f` to each site's neighbourhood using a PairList.
"""
function map_sites!(f::F, out::AbstractVector, nlist::PairList{T,TI}) where {F, T, TI}
    backend = get_array_backend(out)
    nat = nsites(nlist)

    if backend isa CPU
        @threads for i in 1:nat
            j, R = neigs(nlist, i)
            out[i] = f(R)
        end
    else
        # GPU kernel not feasible: f(R) expects dynamically-sized vector,
        # which is incompatible with GPU execution. For GPU workflows,
        # use map_pairs! or map_pairs_d! instead.
        error("map_sites! is not supported on GPU. Use map_pairs! or map_pairs_d! for GPU-accelerated operations.")
    end

    return out
end

"""
    map_pairs!(f, out, nlist::PairList; symmetric=true)

Apply function `f(i, j, R)` to each pair in the neighbour list and accumulate to `out`.

If `symmetric=true`, only processes pairs with i < j and adds f/2 to both sites.
Uses atomic operations when needed for thread safety.
"""
function map_pairs!(f::F, out::AbstractVector, nlist::PairList{T,TI};
                    symmetric::Bool=true) where {F, T, TI}
    backend = get_array_backend(out)
    np = npairs(nlist)

    if np == 0
        return out
    end

    if backend isa CPU
        # For CPU, use thread-local accumulation
        threaded_accumulate!(out, 1:np) do thread_out, n
            i, j = nlist.i[n], nlist.j[n]
            if !symmetric || i < j
                R = _getR(nlist, n)
                val = f(i, j, R)
                if symmetric
                    val_half = val / 2
                    thread_out[i] += val_half
                    thread_out[j] += val_half
                else
                    thread_out[i] += val
                end
            end
        end
    else
        # GPU path: use kernels with atomic accumulation
        if symmetric
            kernel = map_pairs_symmetric_kernel!(backend)
            kernel(f, out, nlist.i, nlist.j, nlist.X, nlist.S, nlist.C; ndrange=np)
        else
            kernel = map_pairs_nonsym_kernel!(backend)
            kernel(f, out, nlist.i, nlist.j, nlist.X, nlist.S, nlist.C; ndrange=np)
        end
        synchronize(backend)
    end

    return out
end

"""
    map_pairs_d!(f, out, nlist::PairList)

Apply anti-symmetric function `f(i, j, R)` to each pair.
Adds `f` to site j and `-f` to site i (for force-like quantities).

For scalar outputs, works on both CPU and GPU.
For vector outputs on GPU, use `map_pairs_d_vec!` with a flat (3, N) matrix.
"""
function map_pairs_d!(f::F, out::AbstractVector{S}, nlist::PairList{T,TI}) where {F, S<:Real, T, TI}
    backend = get_array_backend(out)
    np = npairs(nlist)

    if np == 0
        return out
    end

    if backend isa CPU
        _map_pairs_d_cpu!(f, out, nlist)
    else
        # GPU path: use kernel with atomic accumulation
        kernel = map_pairs_antisym_kernel!(backend)
        kernel(f, out, nlist.i, nlist.j, nlist.X, nlist.S, nlist.C; ndrange=np)
        synchronize(backend)
    end

    return out
end

# Vector version - CPU only for now with SVector output
function map_pairs_d!(f::F, out::AbstractVector{<:SVec}, nlist::PairList{T,TI}) where {F, T, TI}
    backend = get_array_backend(out)
    np = npairs(nlist)

    if np == 0
        return out
    end

    if backend isa CPU
        _map_pairs_d_cpu!(f, out, nlist)
    else
        # GPU with SVector output not directly supported due to atomic limitations
        # Users should use map_pairs_d_vec! with flat (3, N) arrays instead
        error("map_pairs_d! with SVector output not supported on GPU. " *
              "Use map_pairs_d_vec!(f, out_flat, nlist) with a (3, N) matrix instead.")
    end

    return out
end

"""
    map_pairs_d_vec!(f, out_flat::AbstractMatrix, nlist::PairList)

GPU-friendly version of map_pairs_d! for 3D vector outputs.
`out_flat` should be a (3, N) matrix where column i holds the 3D force for atom i.

Example:
```julia
out_flat = CUDA.zeros(Float64, 3, n_atoms)
map_pairs_d_vec!(f, out_flat, nlist)
# Convert back to SVector if needed:
forces = [SVector{3}(out_flat[:, i]) for i in 1:n_atoms]
```
"""
function map_pairs_d_vec!(f::F, out_flat::AbstractMatrix{S}, nlist::PairList{T,TI}) where {F, S<:Real, T, TI}
    backend = get_array_backend(out_flat)
    np = npairs(nlist)
    nat = nsites(nlist)

    @assert size(out_flat, 1) == 3 "out_flat must be (3, N)"
    @assert size(out_flat, 2) == nat "out_flat must have N columns"

    if np == 0
        return out_flat
    end

    if backend isa CPU
        for n in 1:np
            i, j = nlist.i[n], nlist.j[n]
            if i < j
                R = _getR(nlist, n)
                val = f(i, j, R)
                out_flat[1, j] += val[1]
                out_flat[2, j] += val[2]
                out_flat[3, j] += val[3]
                out_flat[1, i] -= val[1]
                out_flat[2, i] -= val[2]
                out_flat[3, i] -= val[3]
            end
        end
    else
        kernel = map_pairs_antisym_vec_kernel!(backend)
        kernel(f, out_flat, nlist.i, nlist.j, nlist.X, nlist.S, nlist.C; ndrange=np)
        synchronize(backend)
    end

    return out_flat
end

# CPU implementation for map_pairs_d! (shared logic)
function _map_pairs_d_cpu!(f::F, out, nlist::PairList) where F
    np = npairs(nlist)
    return threaded_accumulate!(out, 1:np) do thread_out, n
        i, j = nlist.i[n], nlist.j[n]
        if i < j
            R = _getR(nlist, n)
            val = f(i, j, R)
            thread_out[j] += val
            thread_out[i] -= val
        end
    end
end
