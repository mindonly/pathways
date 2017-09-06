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

    # weighted adjacency matrices
    # representing graphs
    #
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

    # use Dijkstra's shortest-path algorithm to compute
    # the minimum energy cost for a graph (adjacency matrix)
    # https://en.wikipedia.org/wiki/Dijkstra's_algorithm
    # [Int]
    #
function mincost(adjmat::Matrix{Int})
    weight = adjmat
    nodect = size(weight)[1]
    source = 1
    target = nodect

    allnodes = Set(1:nodect)
    edges = get_edges(weight)
    distance = Vector{Int}(nodect)
    visited = Set{Int}([source])
    unvisited = setdiff(allnodes, visited)

    distance[source] = 0    # distance from source to source is zero

    while !isempty(unvisited)

            # build Set of candidate paths
        candidates = Set{Tuple{Int, Int}}()
        for i in visited
            for j in unvisited
                if (i, j) in edges
                    push!(candidates, (i, j))
                end
            end
        end

            # if candidate path cost is less than
            # current minimum, update minimum; track node
        min = INFINITY
        node = 0
        for (i, j) in candidates
            if distance[i] + weight[i, j] < min
                min = distance[i] + weight[i, j]
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

    # check if all paths in a path vector
    # span the graph (from source -> target)
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

    # joining tuples in Julia is not built-in
    # https://discourse.julialang.org/t/efficient-tuple-concatenation/5398
    #
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)


    # brute-force connected edge trace,
    # populate results dataframe and display
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

        # collect tuples into Int vectors
        # and finalize
    finalpaths = Vector()
    for path in verifypaths
        a = collect(path)
        push!(finalpaths, a)
    end

        # populate a results dataframe for sorting
    resultsdf = DataFrame(cost = Int[], path = Vector{Int}[])
    for path in finalpaths
        tup = (pathcost(path, graph), unique(path))
        push!(resultsdf, tup)
    end
    resultsdf = sort(resultsdf, cols = :cost, rev = true)

        # display results dataframe
    println("[1] graph traversal paths and costs, \n    brute force connected edge tracing:\n ")
    println(resultsdf, "\n")
end

    # display wrapper for mincost()
    #
function dijkstra_wrapper(adjmat::Matrix{Int})
    graph = adjmat

        # call Dijkstra's mincost()
    println("\n[2] graph minimum cost, \n    Dijkstra's shortest-path:\n ")
    println("\t", mincost(graph), "\n")
end


    # main program
    #
function main()
    println("----------\n GRAPH R: \n----------")
    @time edgetrace(R)
    @time dijkstra_wrapper(R)
    println()
    println("----------\n GRAPH S: \n----------")
    @time edgetrace(S)
    @time dijkstra_wrapper(S)
    println()
    println("----------\n GRAPH T: \n----------")
    @time edgetrace(T)
    @time dijkstra_wrapper(T)
    println()
    println("----------\n GRAPH U: \n----------")
    @time edgetrace(U)
    @time dijkstra_wrapper(U)
    println()
    println("----------\n GRAPH V: \n----------")
    @time edgetrace(V)
    @time dijkstra_wrapper(V)
    println()
    println("----------\n GRAPH W: \n----------")
    @time edgetrace(W)
    @time dijkstra_wrapper(W)

    println("\ntotal runtime: ")
end

@time main()
