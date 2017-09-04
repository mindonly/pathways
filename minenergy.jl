#!/usr/local/bin/env julia

#=
 Rob Sanchez
 CIS 677, F2017
 Project 1
 Calculating Minimum-Energy Reaction Pathways
=#

struct Node
    id::Int
    children::Vector{Int}

    function Node(idx)

        return new(idx, Vector{Int}())
    end
end

struct NodeMap
    nodect::Int
    nodevec::Vector{Node}
    adjmat::Matrix{Int}

    function NodeMap(mat::Matrix{Int})
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

function get_children(node::Node)

    return node.children
end

function get_weight(nmap::NodeMap, idx1::Int, idx2::Int)
    mat = nmap.adjmat
    x = nmap.nodevec[idx1].id
    y = nmap.nodevec[idx2].id

    return get_weight(mat, x, y)
end
get_weight(mat::Matrix{Int}, x::Int, y::Int) = return mat[x, y]

function traverse(nmap::NodeMap, idx::Int)
    if idx != nmap.nodect
        println()
        @show nmap.nodevec[idx].id
        @show nmap.nodevec[idx].children
    end

    energy = 0
    for child in get_children(nmap.nodevec[idx])

        energy = get_weight(nmap, idx, child)
        @show [idx, child], energy

        if child == nmap.nodect
            println("\nPATH COMPLETED!\n")
            break
        end

        energy += traverse(nmap, child)

        @show energy
        println()
    end

    return energy
end

function idxToNode(nmap::NodeMap, idx::Int)

    return nmap.nodevec[idx]
end

# nodemap = NodeMap(R)
# nodemap = NodeMap(S)
# nodemap = NodeMap(U)
nodemap = NodeMap(V)
# nodemap = NodeMap(W)

traverse(nodemap, 1)
