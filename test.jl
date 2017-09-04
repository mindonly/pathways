#!/usr/local/bin/env julia

#=
 Rob Sanchez
 CIS 677, F2017
 Project 1
 Calculating Minimum-Energy Reaction Pathways
=#

struct Node
    id::Int
    next::Vector{Int}

    function Node(idx)
        new(idx, Vector{Int}())
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

function buildmap(matrix)
    dimlen = size(matrix)[1]
    rvec = Vector{Int}()
    for row in 1:dimlen
        for col in 1:dimlen
            if matrix[row, col] != 0
                push!(rvec, row)
            end
        end
    end

    emap = Vector{Node}()
    for row in unique(rvec)
        push!(emap, Node(row))
    end

    for row in 1:dimlen
        for col in 1:dimlen
            if matrix[row, col] != 0
                push!(emap[row].next, col)
            end
        end
    end
    push!(emap, Node(dimlen))


    for n in emap
        @show n.id
        @show n.next
    end
    println()

    return emap
end

function epath(P)
    dimlen = size(P)[1]
    energy = 0
    node = 1
    for row in 1:dimlen
        if node == dimlen
            break
        end
        for col in 1:dimlen
            if node == dimlen
                break
            end
            if P[row, col] != 0
                node = col
                energy += P[row, col]
                @show [row, col], P[row, col], energy, row, node

                if countnz(P[row,:]) > 1
                    P[row, col] = 0
                end

                if row < dimlen
                    row = col
                end
                if col < dimlen
                    col += 1
                end
            end
        end
    end
    println(energy)
    # display(P)
end

# buildmap(R)
# buildmap(S)
# buildmap(U)
# buildmap(V)
# buildmap(W)

adjmat = copy(R)
brs = buildmap(adjmat)
dimlen = size(adjmat)[1]
brslen = length(brs)

energy = 0
for n in brs
    if n.id == brslen
        break
    end
    @show n
    for nxt in n.next
        wt = adjmat[n.id, nxt]
        @show [n.id, nxt], wt
        energy += wt
        n = brs[nxt]
    end
end

println(energy)
