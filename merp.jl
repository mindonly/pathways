#!/usr/local/bin/env julia

#=
 Rob Sanchez
 CIS 677, F2017
 Project 1
 Calculating Minimum-Energy Reaction Pathways
=#

const INFINITY = 9999

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

    wt = get_weight(mat, x, y)

    if wt == 0
        if x == y
            return 0
        else
            return INFINITY
        end
    else
        return wt
    end
end
get_weight(mat::Matrix{Int}, x::Int, y::Int) = return mat[x, y]

function idxToNode(g::Graph, idx::Int)

    return g.nodevec[idx]
end


cost = S
graph = Graph(cost)
n = size(cost)[1]
source = 1
target = n

dist = Vector{Int}(n)
prev = Vector{Node}(n)

for node in get_nodes(graph)
    @show node.id
    @show get_children(node)
    dist[node.id] = get_weight(graph, source, node.id)
    @show source, node.id get_weight(graph, source, node.id)
end

visited = falses(1, n)
visited[source] = true

v = 0
for i in 2:n
    min = INFINITY
    @show min
    for w in 2:n
        if visited[w] == false
            if dist[w] < min
                v = w
                min = dist[w]
            end
        end
    end
    visited[v] = true
    @show v
    for w in 2:n
        if visited[w] == false
            if (min + cost[v, w]) < dist[w]
                dist[w] = min + cost[v, w]
            end
        end
    end
end

println("minimum cost: ", dist[target])
