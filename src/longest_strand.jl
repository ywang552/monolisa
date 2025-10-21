############ Backbone extraction (graph + geometry) ############
# Works with:
#   x = state.x_coords; y = state.y_coords; edges = state.edges
# Assumes edges are undirected pairs (u,v), 1-based node ids.

# --- small helpers ---
const INF = typemax(Float64) / 4



@inline function euclid_len(x::Vector{Float64}, y::Vector{Float64}, u::Int, v::Int)
    dx = x[u] - x[v]; dy = y[u] - y[v]
    return hypot(dx, dy)
end

# Binary heap priority queue (min-heap) for (key, value) = (priority, payload)
mutable struct MinHeap{K,V}
    keys::Vector{K}
    vals::Vector{V}
end
MinHeap{K,V}() where {K,V} = MinHeap{K,V}(K[], V[])
Base.isempty(h::MinHeap) = isempty(h.keys)
@inline function heap_swim!(h::MinHeap, i::Int)
    while i > 1
        p = i >>> 1
        if h.keys[i] < h.keys[p]
            h.keys[i], h.keys[p] = h.keys[p], h.keys[i]
            h.vals[i], h.vals[p] = h.vals[p], h.vals[i]
            i = p
        else
            break
        end
    end
end
@inline function heap_sink!(h::MinHeap, i::Int)
    n = length(h.keys)
    while true
        l = i << 1; r = l + 1; m = i
        if l <= n && h.keys[l] < h.keys[m]; m = l; end
        if r <= n && h.keys[r] < h.keys[m]; m = r; end
        if m == i; break; end
        h.keys[i], h.keys[m] = h.keys[m], h.keys[i]
        h.vals[i], h.vals[m] = h.vals[m], h.vals[i]
        i = m
    end
end
function heappush!(h::MinHeap{K,V}, key::K, val::V) where {K,V}
    push!(h.keys, key); push!(h.vals, val); heap_swim!(h, length(h.keys))
end
function heappop!(h::MinHeap)
    @assert !isempty(h)
    k = h.keys[1]; v = h.vals[1]
    lastk = pop!(h.keys); lastv = pop!(h.vals)
    if !isempty(h)
        h.keys[1] = lastk; h.vals[1] = lastv; heap_sink!(h, 1)
    end
    return k, v
end

# --- graph build ---
"""
    build_graph(x, y, edges)

Returns:
- adj :: Vector{Vector{Int}}  neighbors per node
- w   :: Dict{Tuple{Int,Int}, Float64}  Euclidean length per directed edge (u,v)
"""
function build_graph(x::Vector{Float64}, y::Vector{Float64},
                     edges::Vector{Tuple{Int,Int}}, N::Int)
    adj = [Int[] for _ in 1:N]                 # ensure length == MN
    w   = Dict{Tuple{Int,Int}, Float64}()
    for (u,v) in edges
        # (optional) sanity check:
        # @assert 1 ≤ u ≤ N && 1 ≤ v ≤ N "edge ($(u),$(v)) out of 1..$N"
        l = hypot(x[u]-x[v], y[u]-y[v])
        push!(adj[u], v); push!(adj[v], u)
        w[(u,v)] = l; w[(v,u)] = l
    end
    return adj, w
end


# --- k-core peeling (k=2) on an induced component mask ---
"""
    peel_to_2core!(adj, mask)

Given adjacency and a boolean mask (true=keep), peels degree-1 nodes within the mask.
Returns a tuple (mask_core::Vector{Bool}, removed_stack::Vector{Int}).
"""
function peel_to_2core!(adj::Vector{Vector{Int}}, mask::AbstractVector{Bool})
    n = length(adj)
    deg = zeros(Int, n)
    for u in 1:n
        if mask[u]
            c = 0
            for v in adj[u]; if mask[v]; c += 1; end; end
            deg[u] = c
        end
    end
    Q = Int[]
    for u in 1:n
        if mask[u] && deg[u] <= 1
            push!(Q, u)
        end
    end
    removed = Int[]
    while !isempty(Q)
        u = pop!(Q)
        if !mask[u] || deg[u] > 1; continue; end
        mask[u] = false
        push!(removed, u)
        for v in adj[u]
            if mask[v]
                deg[v] -= 1
                if deg[v] == 1
                    push!(Q, v)
                end
            end
        end
    end
    return mask, removed
end

# --- connected components (over nodes that appear in edges) ---
function components_from_edges(adj::Vector{Vector{Int}})
    n = length(adj)
    seen = falses(n)
    comps = Vector{Vector{Int}}()
    for s in 1:n
        if !seen[s] && !isempty(adj[s])
            q = [s]; seen[s] = true; comp = Int[s]
            while !isempty(q)
                u = pop!(q)
                for v in adj[u]
                    if !seen[v]
                        seen[v] = true; push!(q, v); push!(comp, v)
                    end
                end
            end
            push!(comps, comp)
        end
    end
    return comps
end

# --- PC1 extremes on a set of nodes ---
function pc1_extremes(nodes::Vector{Int}, x::Vector{Float64}, y::Vector{Float64})
    m = length(nodes)
    μx = sum(x[n] for n in nodes) / m
    μy = sum(y[n] for n in nodes) / m
    # covariance of 2D points
    sxx = 0.0; syy = 0.0; sxy = 0.0
    for n in nodes
        dx = x[n]-μx; dy = y[n]-μy
        sxx += dx*dx; syy += dy*dy; sxy += dx*dy
    end
    sxx /= (m-1); syy /= (m-1); sxy /= (m-1)
    C = [sxx sxy; sxy syy]
    evals, evecs = eigen(Symmetric(C))
    # PC1 eigenvector is the one with largest eigenvalue
    idx = (evals[1] > evals[2]) ? 1 : 2
    w = evecs[:, idx]
    # projections
    projmin, amin = +Inf, nodes[1]
    projmax, bmax = -Inf, nodes[1]
    for n in nodes
        proj = w[1]*(x[n]-μx) + w[2]*(y[n]-μy)
        if proj < projmin; projmin = proj; amin = n; end
        if proj > projmax; projmax = proj; bmax = n; end
    end
    return amin, bmax
end

# --- node-based Dijkstra for distances only (for farthest-point pair) ---
function dijkstra_len(adj, w, src::Int)
    n = length(adj)
    dist = fill(INF, n)
    dist[src] = 0.0
    h = MinHeap{Float64,Int}()
    heappush!(h, 0.0, src)
    while !isempty(h)
        d,u = heappop!(h)
        if d > dist[u]; continue; end
        for v in adj[u]
            nd = d + w[(u,v)]
            if nd < dist[v]
                dist[v] = nd
                heappush!(h, nd, v)
            end
        end
    end
    return dist
end

function farthest_pair_in_nodes(adj, w, nodes::Vector{Int})
    s = nodes[1]
    d1 = dijkstra_len(adj, w, s)
    a = argmax(d1[n] for n in nodes)
    a = nodes[a]
    d2 = dijkstra_len(adj, w, a)
    b = argmax(d2[n] for n in nodes)
    b = nodes[b]
    return a, b
end

# --- angle between consecutive directed edges (p->u) then (u->v) ---
@inline function turn_angle(p::Int, u::Int, v::Int, x::Vector{Float64}, y::Vector{Float64})
    ax = x[u]-x[p]; ay = y[u]-y[p]
    bx = x[v]-x[u]; by = y[v]-y[u]
    # normalize
    na = hypot(ax,ay); nb = hypot(bx,by)
    if na == 0.0 || nb == 0.0
        return 0.0
    end
    ca = (ax*bx + ay*by) / (na*nb)
    ca = clamp(ca, -1.0, 1.0)
    return acos(ca)  # radians
end

# --- edge-based (turn-aware) Dijkstra from s to t ---
"""
    turn_aware_path(adj, w, x, y, s, t; λ=0.5)

Returns backbone path (list of node ids) minimizing:
  cost = sum(length(e)) + λ * sum(turn_angle^2)
No turn penalty on the first hop.
"""
function turn_aware_path(adj, w, x, y, s::Int, t::Int; λ::Float64=0.5)
    n = length(adj)

    # State: (prev_node, u) where prev_node==0 means "virtual start at s"
    # Use dictionaries keyed by (prev,u) pairs encoded as Int tuples.
    dist = Dict{Tuple{Int,Int}, Float64}()
    prev = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()  # state -> prior state

    h = MinHeap{Float64,Tuple{Int,Int}}()

    # Initialize from s with no turn penalty
    for v in adj[s]
        st = (s, v)             # came from s to v (prev_node = s, u = v)
        d  = w[(s,v)]
        dist[st] = d
        heappush!(h, d, st)
    end

    best_state = nothing
    best_cost = INF

    while !isempty(h)
        g, (p,u) = heappop!(h)
        if g > get(dist, (p,u), INF); continue; end
        if u == t && g < best_cost
            best_cost = g
            best_state = (p,u)
        end
        # expand to neighbors
        for v in adj[u]
            # avoid immediate backtracking bias if you want (optional):
            # if v == p; continue; end
            ax = x[u]-x[p]; ay = y[u]-y[p]
            bx = x[v]-x[u]; by = y[v]-y[u]
            na = hypot(ax,ay); nb = hypot(bx,by)
            cosθ = (na==0.0 || nb==0.0) ? 1.0 : clamp((ax*bx + ay*by)/(na*nb), -1.0, 1.0)

            # scale so penalty is roughly similar magnitude to θ^2
            λ′ = 0.5 * (π^2) * λ
            cand = g + w[(u,v)] + λ′ * (1 - cosθ)

            st2 = (u, v)
            if cand < get(dist, st2, INF)
                dist[st2] = cand
                prev[st2] = (p,u)
                heappush!(h, cand, st2)
            end
        end
    end

    isnothing(best_state) && return Int[]  # no path

    # Reconstruct nodes from states
    path = Int[t]
    cur = best_state
    while true
        p,u = cur
        pushfirst!(path, u)
        if p == s
            pushfirst!(path, s)
            break
        end
        cur = prev[cur]
    end
    # path may have s duplicated at start; clean consecutive duplicates
    out = Int[]
    last = 0
    for v in path
        if v != last
            push!(out, v)
            last = v
        end
    end
    return out
end

function backbone_for_component(adj, w, x, y, comp::Vector{Int}; λ::Float64=0.5)
    # mask for this component
    n = length(adj)
    mask = falses(n)
    @inbounds for u in comp
        mask[u] = true
    end

    # 2-core peeling
    mask_core, removed = peel_to_2core!(adj, mask)   # accepts BitVector
    core_nodes = [u for u in comp if mask_core[u]]

    # If core empty or tiny, just use farthest-pair on original comp, no turn penalty
    if length(core_nodes) <= 1
        a,b = farthest_pair_in_nodes(adj, w, comp)
        return turn_aware_path(adj, w, x, y, a, b; λ=0.0)
    end

    # Candidate pairs: PC1 extremes + farthest sweeps INSIDE core
    a1, b1 = pc1_extremes(core_nodes, x, y)
    pairs = Tuple{Int,Int}[(a1, b1)]
    append!(pairs, candidate_pairs_farthest(adj, w, core_nodes; K=8))
    pairs = unique(pairs)

    # (Optional) include classic farthest-pair inside core as an extra candidate
    a2, b2 = farthest_pair_in_nodes(adj, w, core_nodes)   # <-- you were missing this
    if (a2, b2) != (a1, b1) && (a2, b2) != (b1, a1)
        push!(pairs, (a2, b2))
    end

    # Evaluate candidates with turn-aware path
    bestP = Int[]; bestScore = -Inf
    for (s, t) in pairs
        P = turn_aware_path(adj, w, x, y, s, t; λ=λ)
        isempty(P) && continue
        # score = geometric length (you can add bend penalty if desired)
        L = 0.0
        @inbounds for i in 1:length(P)-1
            L += w[(P[i], P[i+1])]
        end
        if L > bestScore
            bestScore = L
            bestP = P
        end
    end
    return bestP
end

struct CSR
    rowptr::Vector{Int}   # length = N+1
    colind::Vector{Int}   # length = 2|E| (undirected)
end

"Build undirected CSR from 1-based edges and N (= state.MN)."
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
    return CSR(rowptr, colind)
end

@inline neighbors(csr::CSR, u::Int) = csr.rowptr[u]:(csr.rowptr[u+1]-1)

mutable struct BFSWS
    queue::Vector{Int}
    parent::Vector{Int}   # 0 => no parent
    dist::Vector{Int}     # -1 => unseen
end

init_bfs_ws(N::Int) = BFSWS(Vector{Int}(undef, N), fill(0,N), fill(-1,N))

"Fast reset: only touched nodes are cleared."
@inline function reset_ws!(ws::BFSWS, touched::Vector{Int})
    @inbounds for v in touched
        ws.parent[v] = 0
        ws.dist[v] = -1
    end
    empty!(touched)
end

"Unweighted BFS confined to a component mask (BitVector) if provided."
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

"Double-sweep geodesic confined to a component."
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


"""
    compute_backbone(state; λ=0.5, mode=:both, use_peel=true)

mode:
  :geodesic   -> longest geodesic (no bend penalty, no peel)
  :turnaware  -> straight-ish backbone (2-core + bend penalty)
  :both       -> compute both and return whichever is longer
use_peel controls whether turn-aware step peels to 2-core.
"""
function compute_backbone(state; λ::Float64=0.5, mode=:both, use_peel::Bool=true)
    x = state.x_coords; y = state.y_coords
    N = state.MN
    # Build CSR once
    csr = build_csr_undirected(N, state.edges)

    # Get components (reuse yours; or BFS-on-CSR if you prefer)
    comps = components_all_nodes(N, csr)  # see NOTE below

    # Preallocate one BFS workspace (or one per thread later)
    ws = init_bfs_ws(N)

    result = Dict{Int, Vector{Int}}()
    for (i, comp) in enumerate(comps)
        if length(comp) == 1
            result[i] = comp
            continue
        end

        # 1) Geodesic via double sweep (fast)
        P_geo = geodesic_path_on_comp(csr, ws, comp)

        if mode === :geodesic
            result[i] = isempty(P_geo) ? comp[1:1] : P_geo
            continue
        end

        # 2) Optional: only try turn-aware if geodesic is "short"
        P_turn = Int[]
        if mode !== :geodesic
            if !isempty(P_geo)
                s = P_geo[1]; t = P_geo[end]
                if use_peel
                    P_turn = backbone_for_component_csr(csr, x, y, comp; λ=λ)
                else
                    P_turn = turn_aware_path_csr(csr, x, y, s, t, comp; λ=λ)
                end
            end
        end

        if mode === :turnaware
            result[i] = isempty(P_turn) ? (isempty(P_geo) ? comp[1:1] : P_geo) : P_turn
            continue
        end

        # 3) :both — choose longer by hop-count (unweighted)
        len_geo  = max(length(P_geo)-1, 0)
        len_turn = max(length(P_turn)-1, 0)

        if len_turn > len_geo && !isempty(P_turn)
            result[i] = P_turn
        else
            result[i] = isempty(P_geo) ? comp[1:1] : P_geo
        end
    end
    return result
end


function components_all_nodes(N::Int, csr::CSR)
    seen = falses(N)
    ws = init_bfs_ws(N)
    touched = Int[]
    comps = Vector{Vector{Int}}()
    for s in 1:N
        seen[s] && continue
        # BFS without mask, but capture nodes we touch
        comp = Int[]
        if ws.dist[s] != -1; reset_ws!(ws, touched); end
        far = bfs!(csr, ws, s; touched=touched, compmask=nothing)
        @inbounds for v in touched
            push!(comp, v); seen[v] = true
        end
        push!(comps, comp)
        reset_ws!(ws, touched)
    end
    return comps
end

"""
    candidate_pairs_farthest(adj, w, nodes; K=16)

Return up to K endpoint pairs by repeated farthest-point sweeps
restricted to `nodes`. Great coverage for long/diagonal strands.
"""
function candidate_pairs_farthest(adj, w, nodes::Vector{Int}; K::Int=8)
    pairs = Tuple{Int,Int}[]
    Deg = length.(adj)
    s0 = nodes[argmin(Deg[nodes])]
    seeds = [s0]

    for _ in 1:K
        # farthest sweep from last seed
        d1 = dijkstra_len(adj, w, seeds[end])
        a  = nodes[argmax(d1[n] for n in nodes)]
        d2 = dijkstra_len(adj, w, a)
        b  = nodes[argmax(d2[n] for n in nodes)]
        push!(pairs, (a,b))

        # pick next seed as farthest-from-current-seedset
        # compute distances ONCE per seed, reuse for all n
        Dseeds = [dijkstra_len(adj, w, s) for s in seeds]
        bestn, bestmm = nodes[1], -1.0
        for n in nodes
            # max-min heuristic: distance to nearest seed
            mm = minimum(D[n] for D in Dseeds)
            if mm > bestmm
                bestmm = mm; bestn = n
            end
        end
        push!(seeds, bestn)
        if length(seeds) >= K; break; end
    end
    return unique(pairs)
end

function maximum_spanning_tree_edges_local(adj, w::Dict{Tuple{Int,Int},Float64}, comp::Vector{Int})
    incomp = falses(length(adj))
    @inbounds for u in comp; incomp[u] = true; end

    E = Tuple{Int,Int,Float64}[]
    @inbounds for u in comp
        for v in adj[u]
            # undirected dedup: only keep (u,v) where u < v
            if u < v && incomp[v]
                push!(E, (u, v, w[(u,v)]))
            end
        end
    end
    sort!(E, by = e -> e[3], rev = true)

    # DSU only over comp nodes
    dsu = DisjointSets{Int}()
    for v in comp; push!(dsu, v); end

    T = Vector{Tuple{Int,Int}}()
    for (u,v,_) in E
        if find_root!(dsu, u) != find_root!(dsu, v)
            union!(dsu, u, v)
            push!(T, (u,v))
            if length(T) == length(comp) - 1; break; end
        end
    end
    return T
end

function mst_diameter(adj, w, x, y, comp::Vector{Int})
    Tedges = maximum_spanning_tree_edges_local(adj, w, comp)
    tadj = [Int[] for _ in 1:length(adj)]
    for (u,v) in Tedges
        push!(tadj[u], v); push!(tadj[v], u)
    end
    s = comp[1]
    d1 = dijkstra_len(tadj, w, s)
    a = comp[argmax(d1[n] for n in comp)]
    d2 = dijkstra_len(tadj, w, a)
    b = comp[argmax(d2[n] for n in comp)]
    return turn_aware_path(adj, w, x, y, a, b; λ=0.0)
end


"""
    longest_geodesic(comp)

Return the longest shortest path inside `comp` (weighted diameter).
"""
function longest_geodesic(adj, w, x, y, comp::Vector{Int})
    pairs = candidate_pairs_farthest(adj, w, comp; K=8)
    bestP = Int[]; bestL = -1.0
    for (s,t) in pairs
        P = turn_aware_path(adj, w, x, y, s, t; λ=0.0)  # pure length
        isempty(P) && continue
        L = 0.0
        @inbounds for i in 1:length(P)-1
            L += w[(P[i], P[i+1])]
        end
        if L > bestL
            bestL = L; bestP = P
        end
    end
    return bestP
end
# backbones = compute_backbone(state)