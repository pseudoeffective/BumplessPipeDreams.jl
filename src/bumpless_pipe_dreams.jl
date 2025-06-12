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



############
# droop/drip/drop moves

"""
    can_droop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Check if droop move can be done
"""
function can_droop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if (bpd.mtx[i1,j1],bpd.mtx[i2,j2])  != (2,0)
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2
       if bpd.mtx[i1,j] in [2,3] || bpd.mtx[i2,j] in [2,3]
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2
       if bpd.mtx[i,j1] in [2,3] || bpd.mtx[i,j2] in [2,3]
         return(false)
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.mtx[i,j] in [2,3]
          return(false)
        end
      end
    end
  return(true)
end


# to generate bpds by drooping from flats, need to allow larger droops.
# can have 3 on SW or NE corners of rectangle.

"""
    can_flat_drop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Check if flat drop can be done
"""
function can_flat_drop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

 # check bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # small droop not allowed
    if (i2,j2)==(i1+1,j1+1)
       return(false)
    end

 # check corners
    if (bpd.mtx[i1,j1],bpd.mtx[i2,j2]) != (2,0)
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd.mtx[i1,j] in [2,3,0] || bpd.mtx[i2,j] in [2,3,0]
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd.mtx[i,j1] in [2,3,0] || bpd.mtx[i,j2] in [2,3,0]
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.mtx[i,j] in [2,3,0]
          return(false)
        end
      end
    end
  return(true)

end


# drip is the small droop
"""
    can_drip(bpd::BPD,i1::Int,j1::Int)

Check if a drip (small droop) can be done
"""
function can_drip(bpd::BPD,i1::Int,j1::Int)

 # check corners
    if (bpd.mtx[i1,j1], bpd.mtx[i1+1,j1+1]) != (2,0)
      return(false)
    end

    if !(bpd.mtx[i1,j1+1] in [3,5])
      return(false)
    end

    if !(bpd.mtx[i1+1,j1] in [3,4])
      return(false)
    end

  return(true)

end

# K-theory version
"""
    can_Kdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Check if K-droop move can be done
"""
function can_Kdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if (bpd.mtx[i1,j1],bpd.mtx[i2,j2]) != (2,3)
      return(false)
    end

 # check NE and SW corners
    if !( (bpd.mtx[i1,j2],bpd.mtx[i2,j1]) in [(1,4), (5,1)] )
      return false
    end

 # check N and W borders
    for j=j1+1:j2-1
       if bpd.mtx[i1,j]!=1 && bpd.mtx[i1,j]!=5
         return(false)
       end
    end
    for i=i1+1:i2-1
       if bpd.mtx[i,j1]!=1 && bpd.mtx[i,j1]!=4
         return(false)
       end
    end

 # check S and E borders
    aa=0
    for j=j1+1:j2-1
       if bpd.mtx[i2,j]==3
         return(false)
       end
       if bpd.mtx[i2,j]==2
         aa+=1
         if aa>1 return(false) end
       end
    end

    for i=i1+1:i2-1
       if bpd.mtx[i,j2]==3
         return(false)
       end
       if bpd.mtx[i,j2]==2
         aa+=1
         if aa>1 return(false) end
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.mtx[i,j] in [2,3]
          return(false)
        end
      end
    end
  return(true)
end



# for inverting droop moves


"""
    can_undroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Check if an undroop move can be done
"""
function can_undroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if (bpd.mtx[i1,j1],bpd.mtx[i2,j2])  != (0,3)
      return(false)
    end

# check NE and SW corners
    if !( bpd.mtx[i1,j2] in [2,4] ) || !( bpd.mtx[i2,j1] in [2,5] )
         return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if !( bpd.mtx[i1,j] in [0,4]) || !(bpd.mtx[i2,j] in [1,5])
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if !(bpd.mtx[i,j1] in [0,5]) || !(bpd.mtx[i,j2] in [1,4])
         return(false)
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.mtx[i,j] in [2,3]
          return(false)
        end
      end
    end
  return(true)
end


"""
    can_unKdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Check if an unKdroop move can be done
"""
function can_unKdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

  # check rectangle bounds
  if i2<i1+1 || j2<j1+1
     return false
  end

  # check NW and SE corners
  if (bpd.mtx[i1,j1], bpd.mtx[i2,j2]) != (0,3)
     return false
  end

  # check NE and SW corners
  ne, sw = bpd.mtx[i1,j2], bpd.mtx[i2,j1]
  if !((ne,sw) == (1,2) || (ne,sw) == (2,1))
     return false
  end

  # check N and W borders (at most one 2)
  aa,ii,jj = 0,0,0
  for j=j1+1:j2-1
     if bpd.mtx[i1,j]==3
        return false
     end
     if bpd.mtx[i1,j]==2
        aa+=1
        if aa>1 return false end
     end
  end

  for i=i1+1:i2-1
     if bpd.mtx[i,j1]==3
         return false
     end
     if bpd.mtx[i,j1]==2
        aa+=1
        if aa>1 return false end
     end
  end


  # check S and E borders
  for j=j1+1:j2-1
     if !( bpd.mtx[i2,j] in [1,5] )
         return false
     end
  end
  for i=i1+1:i2-1
     if !( bpd.mtx[i,j2] in [1,4] )
         return false
     end
  end


  # check inside of rectangle
  for i=i1+1:i2-1
    for j=j1+1:j2-1
      if bpd.mtx[i,j] in [2,3]
        return false
      end
    end
  end

  return true
end



"""
    droop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Do droop to BPD on rectangle with NW corner (i1,j1) and SE corner (i2,j2)
"""
function droop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)
 # assumes can_[flat_]droop==true
    
    local bpd2=deepcopy(bpd.mtx)
    # set corners of rectangle
    bpd2[i1,j1]=0
    bpd2[i2,j2]=3
    if bpd2[i1,j2]==5
      bpd2[i1,j2]=2
    elseif bpd2[i1,j2]==3
      bpd2[i1,j2]=4
    end
    if bpd2[i2,j1]==4
      bpd2[i2,j1]=2
    elseif bpd2[i2,j1]==3
      bpd2[i2,j1]=5
    end

    # set west edge
    for i=i1+1:i2-1
       if bpd2[i,j1]==4
          bpd2[i,j1]=0
       elseif bpd2[i,j1]==1
          bpd2[i,j1]=5
       end
    end

    # set north edge
    for j=j1+1:j2-1
       if bpd2[i1,j]==5
          bpd2[i1,j]=0
       elseif bpd2[i1,j]==1
          bpd2[i1,j]=4
       end
    end

    # set east edge
    for i=i1+1:i2-1
       if bpd2[i,j2]==0
          bpd2[i,j2]=4
       elseif bpd2[i,j2]==5
          bpd2[i,j2]=1
       end
    end

    # set south edge
    for j=j1+1:j2-1
       if bpd2[i2,j]==0
          bpd2[i2,j]=5
       elseif bpd2[i2,j]==4
          bpd2[i2,j]=1
       end
    end

    return(BPD(bpd2))
end



"""
    Kdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Do K-droop to BPD on rectangle with NW corner (i1,j1) and SE corner (i2,j2)
"""
function Kdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)
 # assumes can_[flat_]Kdroop==true
    
    local bpd2=deepcopy(bpd.mtx)
    # set corners of rectangle
    bpd2[i1,j1]=0
    bpd2[i2,j2]=3
    if bpd2[i1,j2]==5
      bpd2[i1,j2]=2
    elseif bpd2[i1,j2]==3
      bpd2[i1,j2]=4
    end
    if bpd2[i2,j1]==4
      bpd2[i2,j1]=2
    elseif bpd2[i2,j1]==3
      bpd2[i2,j1]=5
    end

    # set east edge
    ii=i2
    for i=i1+1:i2-1
       if bpd2[i,j2]==0
          bpd2[i,j2]=4
       elseif bpd2[i,j2]==5
          bpd2[i,j2]=1
       elseif bpd2[i,j2]==2
          bpd2[i,j2]=1
          ii=i
       end
    end

    # set south edge
    jj=j2
    for j=j1+1:j2-1
       if bpd2[i2,j]==0
          bpd2[i2,j]=5
       elseif bpd2[i2,j]==4
          bpd2[i2,j]=1
       elseif bpd2[i2,j]==2
          bpd2[i2,j]=1
          jj=j
       end
    end

    # set west edge
    for i=i1+1:min(ii,i2)-1
       if bpd2[i,j1]==4
          bpd2[i,j1]=0
       elseif bpd2[i,j1]==1
          bpd2[i,j1]=5
       end
    end
    # straighten bottom pipe
    if ii<i2
      bpd2[ii,j1]=2
      for j=j1+1:j2-1
        if bpd2[ii,j]==0
          bpd2[ii,j]=5
        elseif bpd2[ii,j]==4
          bpd2[ii,j]=1
        end
      end
    end

    # set north edge
    for j=j1+1:min(jj,j2)-1
       if bpd2[i1,j]==5
          bpd2[i1,j]=0
       elseif bpd2[i1,j]==1
          bpd2[i1,j]=4
       end
    end
    # straighten bottom pipe
    if jj<j2
      bpd2[i1,jj]=2
      for i=i1+1:i2-1
        if bpd2[i,jj]==0
          bpd2[i,jj]=4
        elseif bpd2[i,jj]==5
          bpd2[i,jj]=1
        end
      end
    end

    return(BPD(bpd2))
end

############################
## inverting the droop moves

"""
    undroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Do undroop to BPD on rectangle with NW corner (i1,j1) and SE corner (i2,j2)
"""
function undroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)
    # assumes can_undroop(bpd,i1,j1,i2,j2) == true

    bpd2 = deepcopy(bpd.mtx)

    # restore the four corners
    bpd2[i1,j1] = 2                     # NW
    bpd2[i2,j2] = 0                     # SE

    if bpd2[i1,j2] == 2                 # NE : ╭─ (2) came from ── (5)
        bpd2[i1,j2] = 5
    elseif bpd2[i1,j2] == 4             # NE : │ (4) came from ╯ (3)
        bpd2[i1,j2] = 3
    end

    if bpd2[i2,j1] == 2                 # SW : ╭─ (2) came from │ (4)
        bpd2[i2,j1] = 4
    elseif bpd2[i2,j1] == 5             # SW : ── (5) came from ╯ (3)
        bpd2[i2,j1] = 3
    end

    # west edge
    for i = i1+1 : i2-1
        if bpd2[i,j1] == 0
            bpd2[i,j1] = 4
        elseif bpd2[i,j1] == 5
            bpd2[i,j1] = 1
        end
    end

    # north edge
    for j = j1+1 : j2-1
        if bpd2[i1,j] == 0
            bpd2[i1,j] = 5
        elseif bpd2[i1,j] == 4
            bpd2[i1,j] = 1
        end
    end

    # east edge
    for i = i1+1 : i2-1
        if bpd2[i,j2] == 4
            bpd2[i,j2] = 0
        elseif bpd2[i,j2] == 1
            bpd2[i,j2] = 5
        end
    end

    # south edge
    for j = j1+1 : j2-1
        if bpd2[i2,j] == 5
            bpd2[i2,j] = 0
        elseif bpd2[i2,j] == 1
            bpd2[i2,j] = 4
        end
    end

    return BPD(bpd2)
end



"""
    unKdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)

Do unKdroop to BPD on rectangle with NW corner (i1,j1) and SE corner (i2,j2)
"""
function unKdroop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int)
    # assumes can_unKdroop(bpd,i1,j1,i2,j2) == true

    bpd2 = deepcopy(bpd.mtx)          # work on a copy

    #####################################################################
    # 0.  Locate the special border r-elbows allowed by Kdroop
    #####################################################################
    ii = findfirst(i -> bpd2[i,j1] == 2,  (i1+1):(i2-1))
    if ii === nothing 
       (ii = i2)     # no west‑edge r-elbow
     else ii=ii+i1    # normalize
    end

    jj = findfirst(j -> bpd2[i1,j] == 2,  (j1+1):(j2-1))
    if jj === nothing
        (jj = j2)     # no north‑edge r-elbow
     else jj=jj+j1     # normalize
    end


    #####################################################################
    # 1.  Restore the four corners
    #####################################################################
    bpd2[i1,j1] = 2                         # NW  : 0 → 2
    # NE
    if     bpd2[i1,j2] == 2  bpd2[i1,j2] = 5
    elseif bpd2[i1,j2] == 1  bpd2[i1,j2] = 1
    end
    # SW
    if     bpd2[i2,j1] == 2  bpd2[i2,j1] = 4
    elseif bpd2[i2,j1] == 1  bpd2[i2,j1] = 1
    end
    # SE stays 3  (unchanged by Kdroop)

    #####################################################################
    # 2.  East border  (inverse of: 0→4, 5→1, 2→1)
    #####################################################################

    if ii<i2
      for i = i1+1 : ii-1
          v = bpd2[i,j2]
          if      v == 4                 bpd2[i,j2] = 0
          elseif  v == 1                 bpd2[i,j2] = 5
          end
      end

      bpd2[ii,j2] = 2
    end

    #####################################################################
    # 3.  South border  (inverse of: 0→5, 4→1, 2→1)
    #####################################################################

    if jj<j2
      for j = j1+1 : jj-1
          v = bpd2[i2,j]
          if      v == 5                 bpd2[i2,j] = 0
          elseif  v == 1                 bpd2[i2,j] = 4
          end
      end

      bpd2[i2,jj] = 2
    end


    #####################################################################
    # 4.  West border  (inverse of: 4→0, 1→5)
    #####################################################################
    for i = i1+1 : i2-1
        v = bpd2[i,j1]
        if      v == 0                 bpd2[i,j1] = 4
        elseif  v == 5                 bpd2[i,j1] = 1
    end
    end
    if ii < i2 && ii>i1                      # undo the injected 2
        bpd2[ii,j1] = 4
    end

    #####################################################################
    # 5.  Row ii across the rectangle  (inverse of: 0→5, 4→1)
    #####################################################################
    if ii < i2
        for j = j1+1 : j2-1
            v = bpd2[ii,j]
            if      v == 5             bpd2[ii,j] = 0
            elseif  v == 1             bpd2[ii,j] = 4
            end
        end
    end

    #####################################################################
    # 6.  North border  (inverse of: 5→0, 1→4)
    #####################################################################
    for j = j1+1 : j2-1
        v = bpd2[i1,j]
        if      v == 0                 bpd2[i1,j] = 5
        elseif  v == 4                 bpd2[i1,j] = 1
        end
    end
    if jj < j2 && jj > j1                      # undo the injected 2
        bpd2[i1,jj] = 5
    end

    #####################################################################
    # 7.  Column jj down the rectangle  (inverse of: 0→4, 5→1)
    #####################################################################
    if jj < j2
        for i = i1+1 : i2-1
            v = bpd2[i,jj]
            if      v == 4             bpd2[i,jj] = 0
            elseif  v == 1             bpd2[i,jj] = 5
            end
        end
    end

    return BPD(bpd2)
end


"""
    all_droops(bpd::BPD)

Return vector of all droops of bpd
"""
function all_droops(bpd::BPD)
   local n=size(bpd)[1]

   local dps = BPD[]

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_droop(bpd,i1,j1,i2,j2)
              local bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end



"""
    flat_drops(bpd::BPD)

Return vector of all flat drops of bpd
"""
function flat_drops(bpd::BPD)
   local n=size(bpd.mtx)[1]

   local dps = BPD[]

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_flat_drop(bpd,i1,j1,i2,j2)
              bpd2=makeflat(droop(bpd,i1,j1,i2,j2))
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


"""
    all_drips(bpd::BPD)

Return vector of all drips (small droops) of bpd
"""
function all_drips(bpd::BPD)
   local n=size(bpd.mtx)[1]

   local dps = BPD[]

   for i1=1:n-1
     for j1=1:n-1
            if can_drip(bpd,i1,j1)
              local bpd2=droop(bpd,i1,j1,i1+1,j1+1)
              push!(dps,bpd2)
       end
     end
   end

   return(dps)
end


"""
    all_Kdroops(bpd::BPD)

Return vector of all K-droops of bpd
"""
function all_Kdroops(bpd::BPD)
   local n=size(bpd)[1]

   local dps = Matrix{Int8}[]

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_droop(bpd,i1,j1,i2,j2)
              local bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2.mtx)
            end
            if can_Kdroop(bpd,i1,j1,i2,j2)
              local bpd2=Kdroop(bpd,i1,j1,i2,j2)
              if !(bpd2.mtx in dps)
                push!(dps,bpd2.mtx)
              end
            end
          end
       end
     end
   end

   return(map(BPD,dps))
end





###############
# iterator generating all BPDs for w


struct BPDIterator
    stack::Vector{Tuple{BPD,Vector{BPD}}}
    seen::Set{Matrix}
end


function BPDIterator(bpd::BPD)
    # initialize with the first element
    seen = Set([bpd.mtx])
    droops = all_droops(bpd)
    stack = [(bpd, droops)]
    return BPDIterator(stack,seen)
end

Base.eltype(::Type{BPDIterator}) = BPD

Base.IteratorSize(::Type{<:BPDIterator}) = Base.SizeUnknown()

function Base.iterate(iter::BPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, droops = pop!(iter.stack)

        unseen_droops = filter( b -> !(b.mtx in iter.seen), droops )

        for b in unseen_droops
          push!(iter.seen, b.mtx)  # mark new droop as seen
          push!( iter.stack, (b, all_droops(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # end of iteration
end

########

"""
    all_bpds(w::Vector{Int})

An iterator generating all reduced BPDs for a permutation `w`

## Argument
`w::Vector{Int}`: a permutation

## Returns
`BPDIterator`: an iterator type generating all reduced BPDs for `w`.

## Examples

```julia-repl
# Define the iterator
julia> w = [3,2,5,1,4];

julia> bps = all_bpds(w);

# Run a loop over the iterator

julia> i=0; for b in bps i+=1; end; i
3

# The iterator is exhausted; to reset it, define it again

julia> bps = all_bpds(w);

# Form a vector of all BPDs for w

julia> bpds = collect(bps)
3-element Vector{BPD}:

□ □ ╭─────
□ ╭─┼─────
□ │ │ □ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ ╭─────
□ □ │ ╭───
□ ╭─┼─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ □ ╭───
□ ╭───┼───
□ │ ╭─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─
```
"""
function all_bpds(w::Vector{Int})
    local bpd = Rothe(w)

    iter = BPDIterator(bpd)

    return iter
end


##########
# iterator generating all K-BPDs for w

struct KBPDIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end


function KBPDIterator(bpd::BPD)
    # Initialize with the first element
    seen = Set([bpd.mtx])
    droops = all_Kdroops(bpd)
    stack = [(bpd, droops)]
    return KBPDIterator(stack,seen)
end

Base.eltype(::Type{KBPDIterator}) = BPD

Base.IteratorSize(::Type{<:KBPDIterator}) = Base.SizeUnknown()


function Base.iterate(iter::KBPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, droops = pop!(iter.stack)

        unseen_droops = filter( b -> !(b.mtx in iter.seen), droops )

        for b in unseen_droops
          push!(iter.seen, b.mtx)  # mark new droop as seen
          push!( iter.stack, (b, all_Kdroops(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end



"""
    all_Kbpds(w::Vector{Int})

An iterator generating all BPDs for a permutation `w`, including non-reduced K-theoretic ones

## Argument
`w::Vector{Int}`: a permutation

## Returns
`KBPDIterator`: an iterator type generating all BPDs for `w`.

## Examples
```julia-repl
# Define the iterator
julia> w = [3,2,5,1,4]

julia> bps = all_Kbpds(w);

# Run a loop over the iterator

julia> i=0; for b in bps i+=1; end; i
4

# To reset the iterator, define it again

julia> bps = all_Kbpds(w);

# Form a vector of all BPDs for w

julia> bpds = collect(bps)
4-element Vector{BPD}:

□ □ ╭─────
□ ╭─┼─────
□ │ │ □ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ ╭─────
□ □ │ ╭───
□ ╭─┼─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ □ ╭───
□ □ ╭─┼───
□ ╭─┼─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ □ ╭───
□ ╭───┼───
□ │ ╭─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─ 


```
"""
function all_Kbpds(w::Vector{Int})
    local bpd = Rothe(w)
    iter = KBPDIterator(bpd)

    return iter
end

#######


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
    makeflat(bpd::BPD)

Return the flat BPD in the drift class of bpd
"""
function makeflat(bpd::BPD)
   local n=size(bpd.mtx)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd.mtx[i,j]==0 && bpd.mtx[i-1,j]!=0 && bpd.mtx[i,j-1]!=0 && bpd.mtx[i-1,j-1]==2
         local bpd2=droop(bpd,i-1,j-1,i,j)
         return makeflat(bpd2)
       end
     end
   end

   return(bpd)   

end   
 


#############
# iterator generating all flat BPDs for w

struct FlatBPDIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end

function FlatBPDIterator(bpd::BPD)
    # initialize with the first element
    seen = Set([makeflat(bpd).mtx])
    drops = flat_drops(makeflat(bpd))
    stack = [(makeflat(bpd), drops)]
    return FlatBPDIterator(stack,seen)
end

Base.eltype(::Type{FlatBPDIterator}) = BPD

Base.IteratorSize(::Type{<:FlatBPDIterator}) = Base.SizeUnknown()


function Base.iterate(iter::FlatBPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, drops = pop!(iter.stack)

        unseen_drops = filter( b -> !(makeflat(b).mtx in iter.seen), drops )

        for b in unseen_drops
          b=makeflat(b)
          push!(iter.seen, b.mtx)  # mark new drop as seen
          push!( iter.stack, (b, flat_drops(b)) )
        end

        return( makeflat(current), isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # end of iteration
end


"""
    flat_bpds(w)

An iterator generating all flat reduced BPDs for a permutation `w`

## Argument
`w::Vector{Int}`: a permutation

## Returns
`FlatBPDIterator`: an iterator type generating all flat reduced BPDs for `w`.

## Examples
```julia-repl
# Define the iterator
julia> w = [3,2,5,1,4];

julia> fbps = flat_bpds(w);

# Run a loop over the iterator

julia> i=0; for b in fbps i+=1 end; i
3

# To reset the iterator, define it again

julia> fbps = flat_bpds(w);

# Form a vector of flat BPDs for w

julia> fbpds = collect(fbps)
3-element Vector{BPD}:

□ □ ╭─────
□ ╭─┼─────
□ │ │ □ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ ╭─────
□ □ │ ╭───
□ ╭─┼─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─


□ □ □ ╭───
□ ╭───┼───
□ │ ╭─╯ ╭─
╭─┼─┼───┼─
│ │ │ ╭─┼─
```
"""
function flat_bpds(w)
    local bpd = Rothe(w)
    iter = FlatBPDIterator(bpd)

    return iter
end


function vec_flat_bpds(w)
    local bpd = Rothe(w)
    iter = FlatBPDIterator(bpd)

    return collect(iter)
end



####################
# Words and permutations
####################

function get_sk_from_cross( bpd::BPD, (i,j)::Tuple{Int,Int} )

   if bpd.mtx[i,j] != 1 return nothing end

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




####################
#not used
#=
function can_sharp_drop(bpd,i1,j1,i2,j2)

 # check bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # small droop not allowed
    if (i2,j2)==(i1+1,j1+1)
       return(false)
    end

 # check corners
    if bpd.mtx[i1,j1] != 2 || bpd.mtx[i2,j2] != 0
      return(false)
    end

 # check active pipe
    if i2==i1+1 && bpd.mtx[i2-1,j2-1] == 5
      return(false)
    end
    if j2==j1+1 && bpd.mtx[i2-1,j2-1]==4
      return(false)
    end

 # check NW of destination
    if bpd.mtx[i2-1,j2-1]==0
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd.mtx[i1,j]==2 || bpd.mtx[i2,j]==2 || bpd.mtx[i1,j]==3 || bpd.mtx[i2,j]==3
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd.mtx[i,j1]==2 || bpd.mtx[i,j2]==2 || bpd.mtx[i,j1]==3 || bpd.mtx[i,j2]==3
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.mtx[i,j] == 2 || bpd.mtx[i,j]==3
          return(false)
        end
      end
    end
  return(true)
end



function sharp_drops(bpd::BPD)
# produce all (sharp) drops of bpd
   local n=size(bpd.mtx)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_sharp_drop(bpd,i1,j1,i2,j2)
              bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end



function issharp(bpd::BPD)
# determine if bpd is sharp
   local n=size(bpd.mtx)[1]

   for i=1:n-1
     for j=1:n-1
       if bpd.mtx[i,j]==0 && bpd.mtx[i+1,j]!=0 && bpd.mtx[i,j+1]!=0 && bpd.mtx[i+1,j+1]!=1
          return(false)
       end
     end
   end
   return(true)
end



function makesharp(bpd::BPD)
# returns sharp bpd in drift class of bpd
   local n=size(bpd.mtx)[1]

   for i=1:n-2
     for j=1:n-2
       if bpd.mtx[i,j]==0 && bpd.mtx[i+1,j]!=0 && bpd.mtx[i,j+1]!=0 && bpd.mtx[i+1,j+1]==3
         local bpd2=deepcopy(bpd.mtx)
         bpd2[i,j]=2
         bpd2[i+1,j+1]=0

         if bpd2[i+1,j]==2
           bpd2[i+1,j]=4
         else
           bpd2[i+1,j]=3
         end

         if bpd2[i,j+1]==2
           bpd2[i,j+1]=5
         else
           bpd2[i,j+1]=3
         end

         return makesharp(BPD(bpd2))
       end
     end
   end

   return(bpd)   

end

=#