# Droop moves, etc
# David Anderson, June 2025.


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
    can_flat_drop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int;skew::Bool)

Check if flat drop can be done
"""
function can_flat_drop(bpd::BPD,i1::Int,j1::Int,i2::Int,j2::Int; skew::Bool=false)

 # check that it is not the top pipe
    if skew
      istop(bpd,(i1,j1)) && return(false)
    end


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
    flat_drops(bpd::BPD; skew::Bool)

Return vector of all flat drops of bpd.
If skew==true then don't allow the top pipe to droop.
"""
function flat_drops(bpd::BPD; skew::Bool=false)
   n=size(bpd.mtx)[1]

   dps = BPD[]

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_flat_drop(bpd,i1,j1,i2,j2; skew=skew)
              bpd2=makeflat(droop(bpd,i1,j1,i2,j2); skew=skew)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end

"""
    top_drops(bpd::BPD)

Return vector of all flat drops of bpd obtained by moving only the top pipe.
"""
function top_drops(bpd::BPD)
   n=size(bpd.mtx)[1]

   dps = BPD[]

   for i1=1:n-1
     for j1=1:n-1
       if istop(bpd,(i1,j1))
          for i2=i1+1:n
             for j2=j1+1:n
                 if ( (i2,j2)==(i1+1,j1+1) && can_drip(bpd,i1,j1) ) || can_flat_drop(bpd,i1,j1,i2,j2) 
                    bpd2=makeflat(droop(bpd,i1,j1,i2,j2), skew=true)
                    push!(dps,bpd2)
                 end
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
function flat_bpds(w::Vector{Int})
    local bpd = Rothe(w)
    iter = FlatBPDIterator(bpd)

    return iter
end



#############
# iterator generating all top BPDs for w

struct TopBPDIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end

function TopBPDIterator(bpd::BPD)
    # initialize with the first element
    seen = Set([makeflat(bpd,skew=true).mtx])
    drops = top_drops(makeflat(bpd,skew=true))
    stack = [(makeflat(bpd,skew=true), drops)]
    return TopBPDIterator(stack,seen)
end

Base.eltype(::Type{TopBPDIterator}) = BPD

Base.IteratorSize(::Type{<:TopBPDIterator}) = Base.SizeUnknown()


function Base.iterate(iter::TopBPDIterator, state=nothing)

    while !isempty(iter.stack)
        current, drops = pop!(iter.stack)

        unseen_drops = filter( b -> !(makeflat(b,skew=true).mtx in iter.seen), drops )

        for b in unseen_drops
          b=makeflat(b,skew=true)
          push!(iter.seen, b.mtx)  # mark new drop as seen
          push!( iter.stack, (b, top_drops(b)) )
        end

        return( makeflat(current,skew=true), isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # end of iteration
end

"""
    top_bpds(w)

An iterator generating all top reduced BPDs for a permutation `w`

## Argument
`w::Vector{Int}`: a permutation

## Returns
`TopBPDIterator`: an iterator type generating all top reduced BPDs for `w`.

## Examples
```julia-repl
# Define the iterator
julia> w = [3,1,5,2,4];

julia> tbps = top_bpds(w);

# Run a loop over the iterator

julia> i=0; for b in tbps i+=1 end; i
5

# To reset the iterator, define it again

julia> tbps = top_bpds(w);

# Form a vector of top BPDs for w

julia> collect(tbps)
5-element Vector{BPD}:
                                         
□ □ ╭─────
╭───┼─────
│ □ │ □ ╭─
│ ╭─┼───┼─
│ │ │ ╭─┼─

                
□ □ ╭─────
□ ╭─┼─────
╭─╯ │ □ ╭─
│ ╭─┼───┼─
│ │ │ ╭─┼─

                
□ □ ╭─────
□ □ │ ╭───
╭───┼─╯ ╭─
│ ╭─┼───┼─
│ │ │ ╭─┼─

                
□ □ □ ╭───
□ ╭───┼───
╭─╯ ╭─╯ ╭─
│ ╭─┼───┼─
│ │ │ ╭─┼─

                
□ □ □ ╭───
╭─────┼───
│ □ ╭─╯ ╭─
│ ╭─┼───┼─
│ │ │ ╭─┼─
```
"""
function top_bpds(w::Vector{Int})
    local bpd = Rothe(w)
    iter = TopBPDIterator(bpd)

    return iter
end