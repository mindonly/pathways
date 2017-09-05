#!/usr/local/bin/env julia

#=
 Rob Sanchez
 CIS 677, F2017
 Project 1
 Calculating Minimum-Energy Reaction Pathways
=#

const INFINITY = 9999   # avoid Julia's Inf, which is a Float

mutable struct Node
    id::Int
    visited::Bool
    children::Vector{Int}

    function Node(idx)

        return new(idx, false, Vector{Int}())
    end
end

struct Graph
    nodect::Int
    nodevec::Vector{Node}
    adjmat::Matrix{Int}

    function Graph(mat::Matrix{Int})
        dimlen = size(mat)[1]
        rvec = Vector{Int}()
        for row in 1:dimlen
            for col in 1:dimlen
                if mat[row, col] != 0
                    push!(rvec, row)
                end
            end
        end

        nvec = Vector{Node}()
        for row in unique(rvec)
            push!(nvec, Node(row))
        end

        for row in 1:dimlen
            for col in 1:dimlen
                if mat[row, col] != 0
                    push!(nvec[row].children, col)
                end
            end
        end
        push!(nvec, Node(dimlen))

        # for n in nvec
        #     @show n.id
        #     @show n.children
        # end
        # println()

        return new(dimlen, nvec, mat)
    end
end

R = [0 4 5 0;
     0 0 0 4;
     0 0 0 3;
     0 0 0 0;]

S = [0 4  0  0;
     0 0 10 15;
     0 0  0  3;
     0 0  0  0;]

T = [0 6 0 0  0;
     0 0 2 8 10;
     0 0 0 0 10;
     0 0 0 0  3;
     0 0 0 0  0;]

U = [0 10 0  0  0  0;
     0  0 8 13 24 51;
     0  0 0 14  0  0;
     0  0 0  0  9  0;
     0  0 0  0  0 17;
     0  0 0  0  0  0;]

V = [0 15 14 9  0  0  0;
     0	0  0 0  0 20 37;
     0	5  0 0 17 30  0;
     0	0  0 0 23  0  0;
     0	0  0 0  0  3 20;
     0	0  0 0  0  0 16;
     0	0  0 0  0  0  0;]

W = [0 2 2 2 2 2  0 0  0  0  0  0  0 0;
     0 0 0 0 0 0 10 0  0  0  0  0  0 0;
     0 0 0 0 0 0 12 5  0  0  0  0  0 0;
     0 0 0 0 0 0  0 8 14  0  0  0  0 0;
     0 0 0 0 0 0  0 0  7 11  0  0  0 0;
     0 0 0 0 0 0  0 0  0  2  0  0  0 0;
     0 0 0 0 0 0  0 0  0  0  3  0  0 0;
     0 0 0 0 0 0  0 0  0  0 15  6  0 0;
     0 0 0 0 0 0  0 0  0  0  0 20 13 0;
     0 0 0 0 0 0  0 0  0  0  0  0 18 0;
     0 0 0 0 0 0  0 0  0  0  0  0  0 5;
     0 0 0 0 0 0  0 0  0  0  0  0  0 5;
     0 0 0 0 0 0  0 0  0  0  0  0  0 5;
     0 0 0 0 0 0  0 0  0  0  0  0  0 0;]

function get_nodes(g::Graph)

    return g.nodevec
end

function get_children(node::Node)

    return node.children
end

function get_weight(g::Graph, idx1::Int, idx2::Int)
    mat = g.adjmat
    x = g.nodevec[idx1].id
    y = g.nodevec[idx2].id

    return get_weight(mat, x, y)
end
get_weight(mat::Matrix{Int}, x::Int, y::Int) = return mat[x, y]

function idxToNode(g::Graph, idx::Int)

    return g.nodevec[idx]
end

function get_edges(mat::Matrix{Int})
    dimlen = size(mat)[1]
    edgeset = Set{Tuple{Int, Int}}()
    # edgeset = Set{Vector{Int}}()

    for i = 1:dimlen
        for j = 1:dimlen
            if mat[i, j] != 0
                push!(edgeset, (i, j))
                # push!(edgeset, [i, j])
            end
        end
    end

    return edgeset
end

function mincost(adjmat::Matrix{Int})
    cost = adjmat
    nodect = size(cost)[1]
    source = 1
    target = nodect

    allnodes = Set(1:nodect)
    edges = get_edges(cost)
    distance = Vector{Int}(nodect)
    visited = Set{Int}([source])
    unvisited = setdiff(allnodes, visited)

    distance[source] = 0

    while !isempty(unvisited)
        # Step 1
        candidates = Set{Tuple{Int, Int}}()
        # candidates = Set{Vector{Int}}()
        for i in visited
            for j in unvisited
                if (i, j) in edges
                    push!(candidates, (i, j))
                    # push!(candidates, [i, j])
                    # @show candidates
                end
            end
        end

        # Step 2
        min = INFINITY
        node = 0
        for (i, j) in candidates
            if distance[i] + cost[i, j] < min
                min = distance[i] + cost[i, j]
                node = j
            end
        end

        # Step 3
        distance[node] = min
        push!(visited, node)

        # Step 4
        unvisited = setdiff(allnodes, visited)
    end

    return distance[target]
end

tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

cost = W
edgelist = sort(collect(get_edges(cost)))

source = 1
target = size(cost)[1]

function isconnected(tup1, tup2)
    tuplen = length(tup1)
    if tup1[tuplen] == tup2[1]
        return true
    else
        return false
    end
end

tupvec = Vector()
for item1 in edgelist
    for item2 in edgelist
        if isconnected(item1, item2)
            # println(tuplejoin(item1, item2))
            push!(tupvec, tuplejoin(item1, item2))
        end
    end
end

tupvec2 = Vector()
for item1 in tupvec
    for item2 in edgelist
        if isconnected(item1, item2)
            # println(tuplejoin(item1, item2))
            push!(tupvec2, tuplejoin(item1, item2))
        end
    end
end

tupvec3 = Vector()
for item1 in tupvec2
    for item2 in edgelist
        if isconnected(item1, item2)
            # println(tuplejoin(item1, item2))
            push!(tupvec3, tuplejoin(item1, item2))
        end
    end
end

candidatepaths = vcat(vcat(tupvec, tupvec2), tupvec3)

tuplepaths = Vector()
for item in candidatepaths
    if item[1] == source && item[end] == target
        # println(unique(item))
        push!(tuplepaths, item)
    end
end

finalpaths = Vector()
for path in tuplepaths
    a = collect(path)
    # @show a
    push!(finalpaths, a)
end
# display(finalpaths)

function pathcost(path)
    xvec = Vector()
    yvec = Vector()
    energy = 0
    for (i, elem) in enumerate(path)
        if i % 2 == 0
            push!(yvec, elem)
        else
            push!(xvec, elem)
        end
    end
    for i in zip(xvec, yvec)
        e = collect(i)
        # @show e
        # @show cost[e[1], e[2]]
        energy += cost[e[1], e[2]]
    end

    return energy
end

for path in finalpaths
    println(unique(path), "\t", pathcost(path))
end
println("\nminimum cost: ", mincost(cost))
