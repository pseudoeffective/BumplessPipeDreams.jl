# Tools for ASMs in Julia
# David Anderson, June 2025.




#############
# Conversions to and from ASMs


"""
    bpd2asm(b)

Convert a bumpless pipe dream to an alternating sign matrix

## Argument
`b::BPD`: a bumpless pipe dream

## Returns
`Matrix`: the alternating sign matrix corresponding to `b`.

## Example
```julia-repl
# Generate all reduced BPDs for a permutation
julia> w = [3,2,5,1,4];

julia> bps = all_Kbpds(w);

# Construct a vector of the corresponding ASMs

julia> asms = [];

julia> for b in bps push!(asms, bpd2asm(b)) end;

julia> asms[3]
5×5 Matrix{Int8}:
 0  0  0   1  0
 0  0  1   0  0
 0  1  0  -1  1
 1  0  0   0  0
 0  0  0   1  0

```
"""
function bpd2asm( b::BPD )

  local n=size(b.mtx)[1]

  local a=zeros(Int8,n,n)

  for i=1:n
    for j=1:n
      if b.mtx[i,j]==2
        a[i,j]=1
      elseif b.mtx[i,j]==3
        a[i,j]=-1
      end
    end
  end

  return a
end


"""
    asm2bpd(a)

Convert an alternating sign matrix to a bumpless pipe dream

## Argument
`a::Matrix`: an alternating sign matrix

## Returns
`BPD`: the corresponding BPD

## Example
```julia-repl
# Generate all reduced BPDs for a permutation
julia> w = [3,2,5,1,4];

julia> bpds = collect(all_Kbpds(w));

# Convert a BPD to an ASM

julia> b = bpds[3]
 □ □ □ ╭───
 □ □ ╭─┼───
 □ ╭─┼─╯ ╭─
 ╭─┼─┼───┼─
 │ │ │ ╭─┼─

julia> a = bpd2asm(b)
5×5 Matrix{Int8}:
 0  0  0   1  0
 0  0  1   0  0
 0  1  0  -1  1
 1  0  0   0  0
 0  0  0   1  0

# Convert the ASM back to a BPD

julia> b2 = asm2bpd(a);

julia> b == b2
true

```
"""
function asm2bpd( a::Matrix{<:Integer} )
# improve this, the rules are not local

  local n=size(a)[1]

  local b=Matrix{Int8}(undef,n,n)

  if a[n,1]==1
    b[n,1]=2
  else
    b[n,1]=4
  end

  for j=2:n
    if a[n,j]==1
      b[n,j]=2
    elseif b[n,j-1]==2 || b[n,j-1]==1
      b[n,j]=1
    else
      b[n,j]=4
    end
  end

  for i=n-1:-1:1
    for j=1:n
      if a[i,j]==1
        b[i,j]=2

      elseif a[i,j]==-1
        b[i,j]=3

      elseif a[i,j]==0

        if j==1
          if b[i+1,j]==4
            b[i,j]=4
          else
            b[i,j]=0
          end

        else
          if b[i,j-1]==2 || b[i,j-1]==1 || b[i,j-1]==5
            local bb=5
          else
            bb=0
          end
          if b[i+1,j]==3 || b[i+1,j]==1 || b[i+1,j]==4
            local cc=4
          else
            cc=0
          end
          if bb==5
            if cc==4
              b[i,j]=1
            else
              b[i,j]=5
            end
          else
            b[i,j]=cc
          end
        end
      end
    end
  end

  return BPD(b)
end



"""
    is_asm(mtx::Matrix{<:Integer})

Check if a matrix is a partial ASM, i.e., a NW submatrix of an ASM
"""
function is_asm(mtx::Matrix{<:Integer})
    n,m=size(mtx)
    #check row and column sums
    for i=1:n
        for j=1:m
            if !( sum(mtx[i,1:j]) in [0,1] ) || !( sum(mtx[1:i,j]) in [0,1] ) return false end
        end
    end
    return true
end
