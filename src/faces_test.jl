using Test
# using .Faces  # <- your module from earlier
mutable struct MiniState
    x_coords::Vector{Float64}
    y_coords::Vector{Float64}
    edges::Vector{Tuple{Int,Int}}
end

make_state(xs, ys, edges) = MiniState(collect(xs), collect(ys), collect(edges))

# Point-in-polygon (ray casting), returns true if strictly inside (ignores boundary cases for tests)
function point_in_poly(x::Float64, y::Float64, px::Vector{Float64}, py::Vector{Float64})
    n = length(px); inside = false
    @inbounds for i in 1:n
        j = (i == n) ? 1 : i+1
        xi = px[i]; yi = py[i]
        xj = px[j]; yj = py[j]
        # edge (i->j) crosses horizontal ray at y?
        cond = ((yi > y) != (yj > y)) &&
               (x < (xj - xi) * (y - yi) / (yj - yi + 0.0) + xi)
        if cond; inside = !inside; end
    end
    return inside
end

# --- add: polygon centroid (shoelace-based) ---
function poly_centroid(px::Vector{Float64}, py::Vector{Float64})
    @assert length(px) == length(py) ≥ 3
    A2 = 0.0
    Cx = 0.0
    Cy = 0.0
    n = length(px)
    @inbounds for i in 1:n
        j = (i == n) ? 1 : i+1
        cross = px[i]*py[j] - px[j]*py[i]
        A2 += cross
        Cx += (px[i] + px[j]) * cross
        Cy += (py[i] + py[j]) * cross
    end
    # if area extremely small, fall back to vertex average
    if abs(A2) < 1e-15
        return (sum(px)/n, sum(py)/n)
    end
    A = A2 / 2.0
    Cx /= (3.0*A2)
    Cy /= (3.0*A2)
    return (Cx, Cy)
end

# --- replace the "test point" selection in nesting_parents with centroid ---
function nesting_parents(faces::Vector{Vector{Int}}, x::Vector{Float64}, y::Vector{Float64}, areas::Vector{Float64})
    n = length(faces)
    testpts = Vector{Tuple{Float64,Float64}}(undef, n)
    for i in 1:n
        px = @view x[faces[i]]
        py = @view y[faces[i]]
        testpts[i] = poly_centroid(px, py)  # robust interior point
    end

    parent = fill(0, n)
    for i in 1:n
        (tx, ty) = testpts[i]
        best = 0
        best_area = Inf
        for j in 1:n
            i == j && continue
            px = @view x[faces[j]]; py = @view y[faces[j]]
            if point_in_poly(tx, ty, px, py)   # boundary-free check is fine now
                a = abs(areas[j])
                if a < best_area
                    best_area = a
                    best = j
                end
            end
        end
        parent[i] = best
    end
    return parent
end


# Net region areas = for each top-level boundary loop that is not the global outer boundary,
# subtract areas of its children, add areas of grandchildren, etc. (alternating signs).
function net_region_areas(faces, areas, parent)
    n = length(faces)
    # build children lists
    children = [Int[] for _ in 1:n]
    for i in 1:n
        p = parent[i]
        if p != 0
            push!(children[p], i)
        end
    end
    # recursive alternating-sum
    function alt_sum(i, depth::Int)
        s = (depth % 2 == 0 ? abs(areas[i]) : -abs(areas[i]))
        for c in children[i]
            s += alt_sum(c, depth + 1)
        end
        return s
    end
    # the global outer boundary is the loop with largest |area|
    _, outer_idx = findmax(abs.(areas))
    # compute net area for all top-level loops except the global outer boundary
    nets = Float64[]
    for i in 1:n
        if parent[i] == 0 && i != outer_idx
            push!(nets, alt_sum(i, 0))
        end
    end
    return nets
end

# ---------- TESTS ----------

@testset "Faces: unit square" begin
    xs = [0.0, 1.0, 1.0, 0.0]
    ys = [0.0, 0.0, 1.0, 1.0]
    edges = [(1,2),(2,3),(3,4),(4,1)]
    st = make_state(xs, ys, edges)

    res = Faces.compute_enclosed_faces(st; area_floor=0.0, drop_outer=true,
        normalize_orientation=true, return_abs=true)

    @test length(res.faces) == 1
    @test isapprox(res.areas[1], 1.0; atol=1e-12)
end

@testset "Faces: square + diagonal (two triangles)" begin
    xs = [0.0, 1.0, 1.0, 0.0]
    ys = [0.0, 0.0, 1.0, 1.0]
    edges = [(1,2),(2,3),(3,4),(4,1),(2,4)]  # diagonal splits square
    st = make_state(xs, ys, edges)

    res = Faces.compute_enclosed_faces(st; drop_outer=true,
        normalize_orientation=true, return_abs=true)

    @test length(res.faces) == 2
    total = sum(res.areas)   # two triangles cover the square
    @test isapprox(total, 1.0; atol=1e-12)
end

@testset "Faces: square with inner square hole (net area 3)" begin
    # Outer square 2x2 (area 4), inner 1x1 (area 1), concentric-ish
    xs = [0.0, 2.0, 2.0, 0.0,   0.5, 1.5, 1.5, 0.5]
    ys = [0.0, 0.0, 2.0, 2.0,   0.5, 0.5, 1.5, 1.5]
    # outer loop edges 1-2-3-4-1, inner loop edges 5-6-7-8-5 (no cross edges)
    edges = [(1,2),(2,3),(3,4),(4,1), (5,6),(6,7),(7,8),(8,5)]
    st = make_state(xs, ys, edges)

    # keep signed areas for nesting
    res = Faces.compute_enclosed_faces(st; drop_outer=false,
        normalize_orientation=false, return_abs=false)

    # Build nesting and net region areas, then drop the global outer boundary
    parent = nesting_parents(res.faces, st.x_coords, st.y_coords, res.areas)
    nets = net_region_areas(res.faces, res.areas, parent)
    # There should be exactly one region (the annulus): area 4 - 1 = 3
    @test length(nets) == 1
    @test isapprox(nets[1], 3.0; atol=1e-12)
end

@testset "Orientation sanity: interiors CCW after normalization" begin
    xs = [0.0, 1.0, 1.0, 0.0]
    ys = [0.0, 0.0, 1.0, 1.0]
    edges = [(1,2),(2,3),(3,4),(4,1)]
    st = make_state(xs, ys, edges)

    res = Faces.compute_enclosed_faces(st; drop_outer=false,
        normalize_orientation=true, return_abs=false)

    # largest |area| is outer boundary (negative before normalization or positive after flip)
    _, outer_idx = findmax(abs.(res.areas))
    # All kept interior faces should now be >= 0
    for (i,a) in enumerate(res.areas)
        (i == outer_idx) && continue
        @test a ≥ 0
    end
end
