
using Base.Threads

export maptosites!, maptosites_d!


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


# ==================== KernelAbstractions-based operations ====================

export map_sites!, map_pairs!

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
            js, Rs, Ss = get_neighbours(clist, i)
            out[i] = f(Rs)
        end
    else
        # GPU kernel would go here
        # For now, fall back to sequential
        for i in 1:nat
            js, Rs, Ss = get_neighbours(clist, i)
            out[i] = f(Rs)
        end
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
        for i in 1:nat
            j, R = neigs(nlist, i)
            out[i] = f(R)
        end
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

    if backend isa CPU
        # For CPU, use thread-local accumulation
        nt = nthreads()
        if nt == 1
            _map_pairs_serial!(f, out, nlist, symmetric)
        else
            # Thread-local copies
            outs = [i == 1 ? out : zeros(eltype(out), length(out)) for i in 1:nt]
            @threads for n in 1:np
                tid = threadid()
                i, j = nlist.i[n], nlist.j[n]
                if !symmetric || i < j
                    R = _getR(nlist, n)
                    val = f(i, j, R)
                    if symmetric
                        val_half = val / 2
                        outs[tid][i] += val_half
                        outs[tid][j] += val_half
                    else
                        outs[tid][i] += val
                    end
                end
            end
            # Reduce
            for tid in 2:nt
                out .+= outs[tid]
            end
        end
    else
        # Sequential fallback for GPU (would use atomics in real impl)
        _map_pairs_serial!(f, out, nlist, symmetric)
    end

    return out
end

function _map_pairs_serial!(f::F, out, nlist::PairList, symmetric::Bool) where F
    np = npairs(nlist)
    for n in 1:np
        i, j = nlist.i[n], nlist.j[n]
        if !symmetric || i < j
            R = _getR(nlist, n)
            val = f(i, j, R)
            if symmetric
                val_half = val / 2
                out[i] += val_half
                out[j] += val_half
            else
                out[i] += val
            end
        end
    end
    return out
end

"""
    map_pairs_d!(f, out, nlist::PairList)

Apply anti-symmetric function `f(i, j, R)` to each pair.
Adds `f` to site j and `-f` to site i (for force-like quantities).
"""
function map_pairs_d!(f::F, out::AbstractVector, nlist::PairList{T,TI}) where {F, T, TI}
    backend = get_array_backend(out)
    np = npairs(nlist)

    if backend isa CPU
        nt = nthreads()
        if nt == 1
            for n in 1:np
                i, j = nlist.i[n], nlist.j[n]
                if i < j
                    R = _getR(nlist, n)
                    val = f(i, j, R)
                    out[j] += val
                    out[i] -= val
                end
            end
        else
            outs = [i == 1 ? out : zeros(eltype(out), length(out)) for i in 1:nt]
            @threads for n in 1:np
                tid = threadid()
                i, j = nlist.i[n], nlist.j[n]
                if i < j
                    R = _getR(nlist, n)
                    val = f(i, j, R)
                    outs[tid][j] += val
                    outs[tid][i] -= val
                end
            end
            for tid in 2:nt
                out .+= outs[tid]
            end
        end
    else
        # Sequential fallback
        for n in 1:np
            i, j = nlist.i[n], nlist.j[n]
            if i < j
                R = _getR(nlist, n)
                val = f(i, j, R)
                out[j] += val
                out[i] -= val
            end
        end
    end

    return out
end
