# Tools for generating BPDs in Julia
# David Anderson, June 2025.



#############
# BPD type and constructors


"""
    BPD

A type for bumpless pipe dreams

## Structure

A BPD `b` has one field, `b.mtx`, which is a Matrix{Int8}.  The entries are integers 0-5 which encode the six possible tiles as follows:

   `□ ` <-> `O` <-> `0` (blank)

   `┼─` <-> `+` <-> `1` (cross)

   `╭─` <-> `/` <-> `2` (r-elbow)

   `╯ ` <-> `%` <-> `3` (j-elbow)

   `│ ` <-> `|` <-> `4` (vertical)

   `──` <-> `-` <-> `5` (horizontal)

## Constructor

The function `BPD(m)` takes as its argument `m` either a matrix with entries of type Int8 (with values 0-5) or of type String (with the six possible tiles).


## Example
```julia-repl
julia> mtx = Matrix( [ 0 0 2 ; 0 2 1 ; 2 1 1 ] );

julia> b = BPD( mtx )
 □ □ ╭─
 □ ╭─┼─
 ╭─┼─┼─ 

julia> b == mtx
false

julia> b.mtx == mtx
true

julia> mtx2 = Matrix( [ "O" "O" "/"; "O" "/" "+"; "/" "+" "+" ] );

julia> b2 = BPD( mtx2 )
 □ □ ╭─
 □ ╭─┼─
 ╭─┼─┼─

julia> b2 == b
true

julia> b2.mtx == b.mtx
true
```
"""
struct BPD
    mtx::Matrix{Int8}
end

# Symbol to integer mapping
const SIXVTX_TO_INT = Dict(
    "O" => 0,
    "+" => 1,
    "/" => 2,
    "%" => 3,
    "|" => 4,
    "-" => 5
)

function BPD(matrix::Matrix{String})
    int_matrix = map(x -> SIXVTX_TO_INT[x], matrix)
    return BPD(int_matrix)
end


# convert integers back to symbols for display
# integers 0-5 map to bpd symbols, 6-8 for drifts, 9?
function int_to_symbol(i::Int8)
    symbols = ["\u25A1 ",  # "□ "
        "\u253C\u2500",    # "┼─"
        "\u256D\u2500",    # "╭─"
        "\u256F ",         # "╯ "
        "\u2502 ",         # "│ "
        "\u2500\u2500",    # "──"
        '.', '*', ' ', 'o'
    ]
    return symbols[i+1]
end

function int_to_symbol(t::Tuple)
    return t[1]
end


# add method to Base.show for BPD display
function Base.show(io::IO, bpd::BPD)
    println(io)
    for i in 1:size(bpd.mtx, 1)
        print(" ")
	for j in 1:size(bpd.mtx, 2)
            print(io, int_to_symbol(bpd.mtx[i, j]))
        end
        println(io)
    end
end


# add method to Base.size for BPD
function Base.size(b::BPD)
  size(b.mtx)[1]
end


# overload identity for BPD type
Base.:(==)(b1::BPD, b2::BPD) = b1.mtx == b2.mtx



"""
    Rothe(w::Vector{Int})

Construct the Rothe BPD for a permutation `w``

## Argument
`w::Vector{Int}`: A permutation

## Returns
`BPD`: The Rothe bumpless pipe dream for `w`.

## Example

```julia-repl
# Produce the BPD
julia> w = [3,2,5,1,4];

julia> b = Rothe(w)
 □ □ ╭─────
 □ ╭─┼─────
 □ │ │ □ ╭─
 ╭─┼─┼───┼─
 │ │ │ ╭─┼─


# View the integer matrix which is stored
julia> b.mtx
5×5 Matrix{Int8}:
 0  0  2  5  5
 0  2  1  5  5
 0  4  4  0  2
 2  1  1  5  1
 4  4  4  2  1

```
"""
function Rothe(w::Vector{Int})
        local n=length(w)
        local m=minimum(w)
        local r=Matrix{String}(undef,n,n)

             for j = 1:n
                if j+m-1<w[1]
                   r[1,j]="O"
                elseif j+m-1==w[1]
                   r[1,j]="/"
                else
                   r[1,j]="-"
                end
             end

             for i=2:n
               for j=1:n
                 if j+m-1<w[i]
                   if r[i-1,j]=="/" || r[i-1,j]=="|" || r[i-1,j]=="+"
                      r[i,j]="|"
                   else
                      r[i,j]="O"
                   end
                 elseif j+m-1==w[i]
                   r[i,j]="/"
                 else
                   if r[i-1,j]=="|" || r[i-1,j]=="/" || r[i-1,j]=="+"
                     r[i,j]="+"
                   else
                     r[i,j]="-"
                   end
                 end
               end
             end

        return(BPD(r))
end



"""
    dominant_part(b::BPD)

Extract the partition in NW corner of a BPD
"""
function dominant_part(b::BPD)

  bm = b.mtx

  la = Int[]

  i=1
  while i<size(bm)[1]
    j=1
    lai=0
    while bm[i,j]==0
      lai += 1
      j+=1
    end
    if lai==0
      return la
    else
      push!(la,lai)
      i+=1
    end
  end
  return la    

end


#######
# test properties

# determine if bpd is flat
"""
    isflat(bpd::BPD)

Check if bpd is flat
"""
function isflat(bpd::BPD)
   local n=size(bpd.mtx)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd.mtx[i,j]==0 && bpd.mtx[i-1,j]!=0 && bpd.mtx[i,j-1]!=0 && bpd.mtx[i-1,j-1]==2
          return(false)
       end
     end
   end
   return(true)
end


# returns flat bpd in drift class of bpd
"""
    makeflat(bpd::BPD;skew::Bool)

Return the flat BPD in the drift class of bpd
"""
function makeflat(bpd::BPD; skew::Bool=false)
   n=size(bpd.mtx)[1]

   for i=2:n-1
     for j=2:n-1
       if !skew || !istop(bpd,(i-1,j-1))
         if bpd.mtx[i,j]==0 && bpd.mtx[i-1,j]!=0 && bpd.mtx[i,j-1]!=0 && bpd.mtx[i-1,j-1]==2
           bpd2=droop(bpd,i-1,j-1,i,j)
           return makeflat(bpd2; skew=skew)
         end
        end
     end
   end

   return(bpd)   

end   
 
####################
# Pipe counting
####################
"""
  countpipes( bpd::BPD, (i,j)::Tuple{Int,Int})

Count the number of pipes NW of the pipe at (i,j)
"""
function countpipes( bpd::BPD, (i,j)::Tuple{Int,Int})

  if bpd.mtx[i,j]==0 
    return nothing 
  end

  k=1
  for s in 1:min(i,j)-1
     if bpd.mtx[i-s,j-s] in [2,3,4,5]
        k+=1
     elseif bpd.mtx[i-s,j-s]==1
        k+=2
     end
  end
  return k
end

function istop(bpd::BPD, (i,j)::Tuple{Int,Int})
  return( countpipes(bpd,(i,j))==1 )
end


####################
# Words and permutations
####################


function get_sk_from_cross( bpd::BPD, (i,j)::Tuple{Int,Int} )

   if bpd.mtx[i,j] != 1 return nothing end

   return countpipes(i,j)

end


# extract word from bpd
function bpd2word( bpd::BPD )
  n=size(bpd)
  wd = Int[]

  for s=-n+1:0
      for t=0:s+n-1
        k=get_sk_from_cross(bpd,(n-t,n+s-t))
        k !==nothing && push!(wd, k)
      end
  end
  for s=1:n-1
      for t=0:n-s-1
        k=get_sk_from_cross(bpd,(n-s-t,n-t))
        k !==nothing && push!(wd, k)
      end
  end


  return wd

end


# get permutation associated to a BPD
function bpd2perm( bpd::BPD )

  wrd = bpd2word(bpd)

  return word2perm( wrd )

end


function countboxes( b::BPD )

  ct = count(==(0),b.mtx)
  return ct

end


# decide if a BPD is reduced
function isreduced( bpd::BPD )

  ww = bpd2perm(bpd)

  ct = countboxes(bpd)

  return( len(ww)==ct )

end

####
# copied from permtools

# coxeter length of a permutation (or word)
function len( w::Vector{Int} )
  n = length(w)
  a = 0

  for i in 1:n-1
    for j in i+1:n
      if w[i]>w[j]
        a=a+1
      end
    end
  end

  return a
end


# Compute the Demazure (0-Hecke) product of a word in simple reflections
function word2perm( wrd::Vector{Int} )

  length(wrd)==0 && return [1]

  n = maximum(wrd)+1

  ww = collect(1:n)

  for i=1:length(wrd)
    k=wrd[i]
    if ww[k]<ww[k+1]
     ww = vcat( ww[1:k-1], ww[k+1], ww[k], ww[k+2:n] )
    end
  end

  return ww

end