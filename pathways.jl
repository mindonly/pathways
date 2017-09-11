#!/usr/local/bin/env julia

#=
 Rob Sanchez
 CIS 677, F2017
 Project 1
 Calculating Minimum-Energy Reaction Pathways
=#

using DataFrames

    # avoid Julia's `Inf` (infinity), which is a Float
const INFINITY = 9999

    # digraphs represented by
    # weighted adjacency matrices
    #
G1 = [0 4 5 0;
      0 0 0 4;
      0 0 0 3;
      0 0 0 0;]

G2 = [0 4  0  0;
      0 0 10 15;
      0 0  0  3;
      0 0  0  0;]

G3 = [0 6 0 0  0;
      0 0 2 8 10;
      0 0 0 0 10;
      0 0 0 0  3;
      0 0 0 0  0;]

G4 = [0 10 0  0  0  0;
      0  0 8 13 24 51;
      0  0 0 14  0  0;
      0  0 0  0  9  0;
      0  0 0  0  0 17;
      0  0 0  0  0  0;]

G5 = [0 15 14 9  0  0  0;
      0	0  0 0  0 20 37;
      0	5  0 0 17 30  0;
      0	0  0 0 23  0  0;
      0	0  0 0  0  3 20;
      0	0  0 0  0  0 16;
      0	0  0 0  0  0  0;]

G6 = [0 2 2 2 2 2  0 0  0  0  0  0  0 0;
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

G7 = [0 3 2 1 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0  0;
      0 0 0 0 2 7 0 0 0 0 0 0 0 0 0 0 0 0 0  0;
      0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0  0;
      0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0  0;
      0 0 0 0 0 0 0 2 0 0 3 0 0 0 0 0 0 0 0  0;
      0 0 0 0 0 0 3 2 2 0 0 5 0 0 0 0 0 0 0  0;
      0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0  0;
      0 0 0 0 0 0 0 0 0 0 1 2 0 0 0 0 0 0 0  0;
      0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0  0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 8  0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 0 0  0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 2 1 0 0 0 0  0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2 0 0 0  0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 2 0 10;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 10;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  7;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  5;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  7;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0;]


     # return tuple Set representing graph edges
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

    # check if two graph edges (tuples) are connected
    # [Boolean]
    #
function connected(tup1::NTuple, tup2::NTuple)
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
function pathcost(path::Vector{Int}, adjmat::Matrix{Int})
    cost = adjmat
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

    # link connected edges from two edge lists
    # [Vector{Int}]
    #
function link(list1::Vector, list2::Vector)
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

    # joining tuples in Julia is not built-in
    # https://discourse.julialang.org/t/efficient-tuple-concatenation/5398
    #
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

    # check if all paths in a path vector span the graph
    # (source -> target)
    # [Boolean]
    #
function allspan(pathvec::Vector{Any}, graph::Matrix{Int})
    if isempty(pathvec)
        return false
    end
    source = 1
    target = size(graph)[1]

    for path in pathvec
        if path[1] != source || path[end] != target
            return false
        end
    end

    return true
end


    # brute-force connected edge link & trace
    # [Matrix{Int}]
    #
function edgetrace(adjmat::Matrix{Int})
    graph = adjmat

    source = 1
    target = size(graph)[1]
    edgelist = sort(collect(get_edges(graph)))

        # initial connected edge trace
    rawpathvec = link(edgelist, edgelist)
    candpathvec = copy(rawpathvec)

        # continue path-building tracing connected edges
    while !allspan(rawpathvec, graph)
        rawpathvec = link(rawpathvec, edgelist)
        candpathvec = vcat(candpathvec, rawpathvec)
    end

        # filter for complete paths
    verifypaths = Vector()
    for item in candpathvec
        if item[1] == source && item[end] == target
            push!(verifypaths, item)
        end
    end

        # collect tuples into Int vectors and finalize
    finalpaths = Vector()
    for path in verifypaths
        a = collect(path)
        push!(finalpaths, a)
    end

    return finalpaths
end

    # Dijkstra's algorithm; find shortest path and minimum cost
    # https://en.wikipedia.org/wiki/Dijkstra's_algorithm
    # [Tuple{Int, Vector{Int}}]
    #
function dijkstra(adjmat::Matrix{Int})
    weight = adjmat
    nodect = size(weight)[1]
    source = 1
    target = nodect

    allnodes = Set(1:nodect)
    edges = get_edges(weight)
    distance = Vector{Int}(nodect)
    previous = Vector{Int}(nodect)
    visited = Set{Int}([source])
    unvisited = setdiff(allnodes, visited)

    distance[source] = 0    # distance from source to source is zero

    while !isempty(unvisited)

            # build Set of neighbor paths
        neighbors = Set{Tuple{Int, Int}}()
        for i in visited
            for j in unvisited
                if (i, j) in edges
                    push!(neighbors, (i, j))
                end
            end
        end

            # if neighbor path cost is less than current
            # minimum, update minimum; track node & prev
        min = INFINITY
        node = 0
        prev = 0
        for (i, j) in neighbors
            if distance[i] + weight[i, j] < min
                min = distance[i] + weight[i, j]
                node = j
                prev = i
            end
        end

            # update distance to node from source
            # track previous node
            # update visited Set
        distance[node] = min
        previous[node] = prev
        push!(visited, node)

            # update unvisited Set by set difference
        unvisited = setdiff(allnodes, visited)
    end

        # reverse-iterate shortest path
    shortpath = Vector{Int}()
    push!(shortpath, target)
    v = target
    while previous[v] != source
        push!(shortpath, previous[v])
        v = previous[v]
    end
    push!(shortpath, source)

    return (distance[target], reverse(shortpath))
end

    # display wrapper for edgetrace()
    #
function edgetrace_wrapper(adjmat::Matrix{Int})
    println("[1] graph traversal costs and paths,\n    brute force connected edge tracing:\n")

    graph = adjmat
    @time pathv = edgetrace(graph)
    println()

        # populate a results dataframe for sorting
    resultsdf = DataFrame(cost = Int[], path = Vector{Int}[])
    for path in pathv
        tup = (pathcost(path, graph), unique(path))
        push!(resultsdf, tup)
    end
    resultsdf = sort(resultsdf, cols = :cost, rev = true)

        # display results dataframe
    println(resultsdf)
end

    # display wrapper for dijkstra()
    #
function dijkstra_wrapper(adjmat::Matrix{Int})
    graph = adjmat

        # call dijkstra()
    println("\n[2] graph minimum cost and optimal path, \n    Dijkstra's algorithm:\n ")
    @time println("\t", dijkstra(graph), "\n")
end


    # main program
    #
function main()

        # call edgetrace() and dijkstra() on each graph
    for (i, graph) in enumerate([G1, G2, G3, G4, G5, G6, G7])
        println("----------\n GRAPH $i: \n----------")
        edgetrace_wrapper(graph)
        dijkstra_wrapper(graph)
        println()
    end

    println("\n----------------\n total runtime: \n----------------")
end

@time main()
