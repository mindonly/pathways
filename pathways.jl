#!/usr/local/bin/env julia

#=
 Rob Sanchez
 CIS 677, F2017
 Project 1
 Calculating Minimum-Energy Reaction Pathways
=#

const INFINITY = 9999   # avoid Julia's `Inf` (infinity), which is a Float

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


     # return tuple Set representing graph
     # edges (adjacency matrix)
     # [Set{Tuple{Int, Int}}]
     #
function get_edges(mat::Matrix{Int})
    dimlen = size(mat)[1]
    edgeset = Set{Tuple{Int, Int}}()

    for i = 1:dimlen
        for j = 1:dimlen
            if mat[i, j] != 0
                push!(edgeset, (i, j))
            end
        end
    end

    return edgeset
end

    # check if two graph edges (tuples)
    # are connected
    # [Boolean]
    #
function connected(tup1, tup2)
    tuplen = length(tup1)

    if tup1[tuplen] == tup2[1]
        return true
    else
        return false
    end
end

    # return energy cost for a single graph path
    # [Int]
    #
function pathcost(path::Vector{Int})
    xvec = Vector{Int}()
    yvec = Vector{Int}()
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
        energy += cost[e[1], e[2]]
    end

    return energy
end

    # use Dijkstra's shortest-path algorithm to compute
    # the minimum energy cost for a graph (adjacency matrix)
    # https://en.wikipedia.org/wiki/Dijkstra's_algorithm
    # [Int]
    #
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

            # build Set of candidate paths
        candidates = Set{Tuple{Int, Int}}()
        for i in visited
            for j in unvisited
                if (i, j) in edges
                    push!(candidates, (i, j))
                    # @show candidates
                end
            end
        end

            # if candidate path cost is less than
            # current minimum, update
        min = INFINITY
        node = 0
        for (i, j) in candidates
            if distance[i] + cost[i, j] < min
                min = distance[i] + cost[i, j]
                node = j
            end
        end

            # update distance
            # to node from source
        distance[node] = min
        push!(visited, node)

            # update unvisited set by set difference
        unvisited = setdiff(allnodes, visited)
    end

    return distance[target]
end

    # joining Julia tuples is not built-in
    # https://discourse.julialang.org/t/efficient-tuple-concatenation/5398
    #
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)






# MAIN PROGRAM #

cost = W
edgelist = sort(collect(get_edges(cost)))

source = 1
target = size(cost)[1]



function trailblaze(list1, list2)
    tempvec = Vector()
    for item1 in list1
        for item2 in list2
            if connected(item1, item2)
                push!(tempvec, tuplejoin(item1, item2))
            end
        end
    end

    return tempvec
end

function allspan(pathvec)
    if isempty(pathvec)
        return false
    end
    for path in pathvec
        if path[1] != source || path[end] != target
            return false
        end
    end

    return true
end


tupvec = Vector()
for item1 in edgelist
    for item2 in edgelist
        if connected(item1, item2)
            push!(tupvec, tuplejoin(item1, item2))
        end
    end
end

tupvec2 = Vector()
for item1 in tupvec
    for item2 in edgelist
        if connected(item1, item2)
            push!(tupvec2, tuplejoin(item1, item2))
        end
    end
end

tupvec3 = Vector()
for item1 in tupvec2
    for item2 in edgelist
        if connected(item1, item2)
            push!(tupvec3, tuplejoin(item1, item2))
        end
    end
end

candidatepaths = vcat(vcat(tupvec, tupvec2), tupvec3)


tuplepaths = Vector()
for item in candidatepaths
    if item[1] == source && item[end] == target
        push!(tuplepaths, item)
    end
end

finalpaths = Vector()
for path in tuplepaths
    a = collect(path)
    push!(finalpaths, a)
end



for path in finalpaths
    println(unique(path), "\t", pathcost(path))
end
println("\nminimum cost: ", mincost(cost))


candpathvec = trailblaze(edgelist, edgelist)
finalpathvec = copy(candpathvec)
display(candpathvec)

while allspan(candpathvec) == false
    candpathvec = trailblaze(candpathvec, edgelist)
    finalpathvec = vcat(finalpathvec, candpathvec)
end


tuplepaths = Vector()
for item in finalpathvec
    if item[1] == source && item[end] == target
        push!(tuplepaths, item)
    end
end

finalpaths = Vector()
for path in tuplepaths
    a = collect(path)
    push!(finalpaths, a)
end



for path in finalpaths
    println(unique(path), "\t", pathcost(path))
end
println("\nminimum cost: ", mincost(cost))
