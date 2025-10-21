
# compute_backbone_fast.jl
# Self-contained, allocation-lean backbone extraction for very large graphs (~1e6 nodes).
# Dependencies: none (pure Julia). Drop into your project and `include("compute_backbone_fast.jl")`.
#
# Assumes `state` has:
#   - state.x_coords::Vector{Float64}
#   - state.y_coords::Vector{Float64}
#   - state.edges::Vector{Tuple{Int,Int}}   (1-based undirected edges)
#   - state.MN::Int                         (number of nodes)
#
# API:
#   compute_backbone(state; λ=0.5, mode=:both, use_peel=true) -> Dict{Int,Vector{Int}}
#
# Modes:
#   :geodesic  - BFS double-sweep geodesic backbone (very fast)
#   :turnaware - single turn-aware path between geodesic endpoints
#   :both      - compute geodesic; try turn-aware once; return the longer path

############################
# Minimal binary Min-Heap  #
############################

struct MinHeap{T}
    data::Vector{T}
end
MinHeap{T}() where {T} = MinHeap{T}(T[])
Base.isempty(h::MinHeap) = isempty(h.data)

@inline function sift_up!(v)
    i = length(v)
    while i > 1
        p = i >>> 1
        if isless(v[i], v[p])
            v[i], v[p] = v[p], v[i]
            i = p
        else
            break
        end
    end
end

@inline function sift_down!(v, i::Int)
    n = length(v)
    while true
        l = i << 1
        r = l + 1
        if l > n; break; end
        m = (r <= n && isless(v[r], v[l])) ? r : l
        if isless(v[m], v[i])
            v[i], v[m] = v[m], v[i]
            i = m
        else
            break
        end
    end
end

@inline function heappush!(h::MinHeap{T}, x::T) where {T}
    push!(h.data, x)
    sift_up!(h.data)
    return nothing
end

@inline function heappop!(h::MinHeap{T}) where {T}
    @assert !isempty(h)
    v = h.data
    top = v[1]
    v[1] = v[end]
    pop!(v)
    if !isempty(v)
        sift_down!(v, 1)
    end
    return top
end

############################
# CSR graph representation #
############################

struct CSR
    rowptr::Vector{Int}   # length = N+1
    colind::Vector{Int}   # length = 2|E| (undirected stored as both directions)
    src::Vector{Int}      # length = 2|E| ; src[k] = u for edge u -> colind[k]
end

"Build undirected CSR with per-edge source indices."
function build_csr_undirected(N::Int, edges::Vector{Tuple{Int,Int}})
    deg = zeros(Int, N)
    @inbounds for (u,v) in edges
        deg[u] += 1; deg[v] += 1
    end
    rowptr = Vector{Int}(undef, N+1)
    rowptr[1] = 1
    @inbounds for i in 1:N
        rowptr[i+1] = rowptr[i] + deg[i]
    end
    colind = Vector{Int}(undef, rowptr[end]-1)
    fill!(deg, 0) # reuse as cursors
    @inbounds for (u,v) in edges
        pu = rowptr[u] + deg[u]; colind[pu] = v; deg[u] += 1
        pv = rowptr[v] + deg[v]; colind[pv] = u; deg[v] += 1
    end
    src = Vector{Int}(undef, length(colind))
    @inbounds for u in 1:N
        for k in rowptr[u]:(rowptr[u+1]-1)
            src[k] = u
        end
    end
    return CSR(rowptr, colind, src)
end

@inline neighbors(csr::CSR, u::Int) = csr.rowptr[u]:(csr.rowptr[u+1]-1)

###########################
# BFS workspace + helpers #
###########################

const INF = typemax(Float64) / 4

mutable struct BFSWS
    queue::Vector{Int}
    parent::Vector{Int}   # 0 => no parent
    dist::Vector{Int}     # -1 => unseen
end

init_bfs_ws(N::Int) = BFSWS(Vector{Int}(undef, N), fill(0, N), fill(-1, N))

"Reset only touched nodes to avoid O(N) clears."
@inline function reset_ws!(ws::BFSWS, touched::Vector{Int})
    @inbounds for v in touched
        ws.parent[v] = 0
        ws.dist[v] = -1
    end
    empty!(touched)
end

"Unweighted BFS confined by an optional component mask."
function bfs!(csr::CSR, ws::BFSWS, src::Int; touched::Vector{Int}, compmask::Union{Nothing,BitVector}=nothing)
    q = ws.queue; parent = ws.parent; dist = ws.dist
    head = 1; tail = 0

    dist[src] = 0
    push!(touched, src)
    tail += 1; q[tail] = src

    far = src
    @inbounds while head <= tail
        u = q[head]; head += 1
        du = dist[u]
        for k in neighbors(csr, u)
            v = csr.colind[k]
            if compmask !== nothing && !compmask[v]; continue; end
            if dist[v] == -1
                dist[v] = du + 1
                parent[v] = u
                push!(touched, v)
                tail += 1; q[tail] = v
                if dist[v] > dist[far]; far = v; end
            end
        end
    end
    return far
end

"Parent-based path reconstruction."
function reconstruct_path(parent::Vector{Int}, s::Int, t::Int)
    if s == 0 || t == 0; return Int[]; end
    P = Int[]
    v = t
    while v != 0
        push!(P, v)
        v == s && break
        v = parent[v]
    end
    reverse!(P)
    return P
end

"Double-sweep geodesic path confined to a component."
function geodesic_path_on_comp(csr::CSR, ws::BFSWS, comp::Vector{Int})
    if isempty(comp); return Int[]; end
    compmask = falses(length(ws.parent))
    @inbounds for v in comp; compmask[v] = true; end
    touched = Int[]
    a = bfs!(csr, ws, comp[1]; touched=touched, compmask=compmask)
    reset_ws!(ws, touched)
    b = bfs!(csr, ws, a; touched=touched, compmask=compmask)
    P = reconstruct_path(ws.parent, a, b)
    reset_ws!(ws, touched)
    return P
end

###########################
# Connected components    #
###########################

"Connected components via BFS on CSR; returns Vector of node lists."
function components_all_nodes(N::Int, csr::CSR)
    seen = falses(N)
    comps = Vector{Vector{Int}}()
    q = Vector{Int}(undef, N)
    for s in 1:N
        if seen[s]; continue; end
        head = 1; tail = 0
        tail += 1; q[tail] = s; seen[s] = true
        comp = Int[]
        @inbounds while head <= tail
            u = q[head]; head += 1
            push!(comp, u)
            for k in neighbors(csr, u)
                v = csr.colind[k]
                if !seen[v]
                    seen[v] = true
                    tail += 1; q[tail] = v
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end

###########################
# Turn-aware path (CSR)   #
###########################

struct EdgeRec
    cost::Float64
    edge::Int  # k index into csr.colind/src
end
Base.isless(a::EdgeRec, b::EdgeRec) = a.cost < b.cost

"Turn-aware shortest path from s to t over CSR using hop cost + bend penalty."
function turn_aware_path_csr(csr::CSR, x::Vector{Float64}, y::Vector{Float64},
                             s::Int, t::Int; λ::Float64=0.5)
    M = length(csr.colind)
    dist_e   = fill(INF, M)   # cost per directed edge
    parent_e = fill(0,   M)   # previous edge id

    h = MinHeap{EdgeRec}()

    # scale penalty similar to original: λ′ = 0.5 * π^2 * λ
    λ′ = 0.5 * (pi^2) * λ

    # Seed with edges out of s (no bend penalty on first hop)
    @inbounds for k in neighbors(csr, s)
        dist_e[k] = 1.0                     # hop cost
        parent_e[k] = 0
        heappush!(h, EdgeRec(dist_e[k], k))
    end

    best_edge = 0
    best_cost = INF

    @inbounds while !isempty(h)
        rec = heappop!(h)
        k = rec.edge
        g = rec.cost
        if g > dist_e[k]; continue; end  # stale

        u = csr.src[k]
        v = csr.colind[k]

        if v == t && g < best_cost
            best_cost = g
            best_edge = k
            # optional: break  # keeps optimal if you skip stale pops
        end

        for knext in neighbors(csr, v)
            w = csr.colind[knext]

            # turn penalty between (u->v) and (v->w)
            ax = x[v]-x[u]; ay = y[v]-y[u]
            bx = x[w]-x[v]; by = y[w]-y[v]
            na = hypot(ax, ay); nb = hypot(bx, by)
            cosθ = (na==0.0 || nb==0.0) ? 1.0 : clamp((ax*bx + ay*by)/(na*nb), -1.0, 1.0)

            cand = g + 1.0 + λ′*(1 - cosθ)  # hop + bend

            if cand < dist_e[knext]
                dist_e[knext] = cand
                parent_e[knext] = k
                heappush!(h, EdgeRec(cand, knext))
            end
        end
    end

    best_edge == 0 && return Int[]

    # reconstruct nodes from edge chain
    edges = Int[]
    e = best_edge
    while e != 0
        push!(edges, e)
        e = parent_e[e]
    end
    out = Vector{Int}(undef, length(edges)+1)
    out[end] = t
    @inbounds for i in 1:length(edges)
        e = edges[i]
        out[end - i] = csr.src[e]
    end

    # clean adjacent duplicates (rare)
    cleaned = Int[]
    last = 0
    @inbounds for v in out
        if v != last; push!(cleaned, v); last = v; end
    end
    return cleaned
end

###########################
# 2-core peel on a comp   #
###########################

"Return mask (BitVector) of component 2-core (true=kept)."
function peel_to_2core_csr!(csr::CSR, comp::Vector{Int})
    N = length(csr.rowptr) - 1
    incomp = falses(N)
    @inbounds for u in comp; incomp[u] = true; end
    deg = zeros(Int, N)
    @inbounds for u in comp
        c = 0
        for k in neighbors(csr, u)
            v = csr.colind[k]
            c += incomp[v] ? 1 : 0
        end
        deg[u] = c
    end
    Q = Int[]
    @inbounds for u in comp
        if deg[u] <= 1; push!(Q, u); end
    end
    @inbounds while !isempty(Q)
        u = pop!(Q)
        if !incomp[u] || deg[u] > 1; continue; end
        incomp[u] = false
        for k in neighbors(csr, u)
            v = csr.colind[k]
            if incomp[v]
                deg[v] -= 1
                if deg[v] == 1; push!(Q, v); end
            end
        end
    end
    return incomp
end

"Backbone on component: geodesic endpoints -> (optional) 2-core -> turn-aware."
function backbone_for_component_csr(csr::CSR, x::Vector{Float64}, y::Vector{Float64},
                                    comp::Vector{Int}; λ::Float64=0.5, use_peel::Bool=true)
    ws = init_bfs_ws(length(x))
    Pgeo = geodesic_path_on_comp(csr, ws, comp)
    isempty(Pgeo) && return Int[]
    s, t = Pgeo[1], Pgeo[end]

    if !use_peel
        return turn_aware_path_csr(csr, x, y, s, t; λ=λ)
    end

    mask = peel_to_2core_csr!(csr, comp)
    # If core is tiny or s/t not in core, use original endpoints
    core_nodes = Vector{Int}()
    sizehint!(core_nodes, length(comp))
    @inbounds for u in comp
        if mask[u]; push!(core_nodes, u); end
    end
    if length(core_nodes) <= 1 || !(mask[s] && mask[t])
        return turn_aware_path_csr(csr, x, y, s, t; λ=λ)
    end

    # Re-pick endpoints inside core using geodesic on the core
    Pcore = geodesic_path_on_comp(csr, ws, core_nodes)
    isempty(Pcore) && return turn_aware_path_csr(csr, x, y, s, t; λ=λ)
    a, b = Pcore[1], Pcore[end]
    return turn_aware_path_csr(csr, x, y, a, b; λ=λ)
end

###########################
# Top-level entry point   #
###########################

"""
compute_backbone(state; λ=0.5, mode=:both, use_peel=true)

mode:
  :geodesic   -> geodesic (double-sweep BFS)
  :turnaware  -> single turn-aware run (2-core optional)
  :both       -> compute geodesic, try turn-aware once, return longer
"""
function compute_backbone(state; λ::Float64=0.5, mode=:both, use_peel::Bool=true)
    x = state.x_coords::Vector{Float64}
    y = state.y_coords::Vector{Float64}
    N = state.MN::Int
    edges = state.edges::Vector{Tuple{Int,Int}}

    csr = build_csr_undirected(N, edges)
    comps = components_all_nodes(N, csr)

    ws = init_bfs_ws(N)

    result = Dict{Int, Vector{Int}}()
    for (i, comp) in enumerate(comps)
        if length(comp) == 1
            result[i] = comp
            continue
        end

        # 1) Fast geodesic
        P_geo = geodesic_path_on_comp(csr, ws, comp)

        if mode === :geodesic
            result[i] = isempty(P_geo) ? comp[1:1] : P_geo
            continue
        end

        # 2) Single turn-aware attempt (no K-pair scans)
        P_turn = Int[]
        if !isempty(P_geo)
            s, t = P_geo[1], P_geo[end]
            if use_peel
                P_turn = backbone_for_component_csr(csr, x, y, comp; λ=λ, use_peel=true)
            else
                P_turn = turn_aware_path_csr(csr, x, y, s, t; λ=λ)
            end
        end

        if mode === :turnaware
            result[i] = isempty(P_turn) ? (isempty(P_geo) ? comp[1:1] : P_geo) : P_turn
            continue
        end

        # 3) :both — choose longer by hop-count
        len_geo  = max(length(P_geo)-1, 0)
        len_turn = max(length(P_turn)-1, 0)
        result[i] = (len_turn > len_geo && !isempty(P_turn)) ? P_turn :
                    (isempty(P_geo) ? comp[1:1] : P_geo)
    end
    return result
end

# End of file.
