
# This file "ospall" contains the function: ospall().
# Date: 2018-03-30
# Version: 1.01
# Author: Jaap de Gruijter

# CONTENT:
# 1. Systematic sampling from the grid.
# 2. Ospats application to the sample.
# 3. Allocation of non-sample points to sample-strata.
# 4. Sample sizes, total and in strata.
# 5. Stratified random sampling.
# 6. Output of final results.

# For process monitoring uncomment println(...) lines

function ospall()
println("------------- START FUNCTION OSPALL -------------------")

##### SECTION 1. SYSTEMATIC SAMPLING FROM THE GRID WITH RANDOM
##### STARTING POINT
inrow = 1:in
start = rand(rng, inrow)[1]
sample = start:in:N
n = length(sample)
println("sample size = ", n)
x_s = x[sample]
y_s = y[sample]
z_pred_s = z_pred[sample]
s2_s = s2[sample]

##### SECTION 2. APPLICATION OF OSPATS METHOD TO THE SAMPLE
##### Subsection 2.1 Calculate n x n matrix of generalized distances
# println("----- Calculating matrix of generalized distances --")
d2 = zeros(n,n)
for i = 1:(n-1)
  for j = (i+1):n
    d2[i,j] = ((z_pred_s[i] - z_pred_s[j])^2)/R2 + (s2_s[i] + s2_s[j]).*(1 - exp(-3*(sqrt((x_s[i] - x_s[j])^2 + (y_s[i] - y_s[j])^2)/range)))
  end
end
d2 = d2 + d2'

TOTd2 = sum(d2)/2
ObarH1 = sqrt(TOTd2)/n
# println("ObarH1 = ", ObarH1)
strat0 = Array{Int64,1}(n)
stratcy = Array{Int64,1}(n)
stratbest = Array{Int64,1}(n)
global stratbest, nhbest, Hbest, ObarFinal, Obar_grid, nbest, nh_l, d2, n, Nh_sample, Nh, nh

#######################################   Start optimization of H

for H = H_max : -1 : H_min
println("_______________________________ Number of strata : ",H)

cbObj = Array{Float64,1}(H)
cbObj = zeros(Float64,H,1)
TotTransf = 0

##### Subsection 2.3. Initial stratification
  missing = n - H*floor(Int64, n/H)
  A = collect(1:H)
  B = vcat(A,A)
  repeat = n/H -2
  for rep = 1:repeat
    B = vcat(B,A)
  end
  fillup = collect(1:missing)
  B = vcat(B,fillup)

  v = collect(1:n)
  w = v[randperm(rng,n)]
  for i = 1:n
    strat0[w[i]] = B[i]
  end

##### Subsection 2.4. Contributions from sample-strata to O
Sd2 = Array{Float64,1}(H)
Sd2 = zeros(Float64,H,1)

for strat = 1:H
  Sd2[strat] = 0
  for i = 1:(n-1)
    if strat0[i] == strat
      for j = (i+1):n
        if strat0[j] == strat
          Sd2[strat] = Sd2[strat] + d2[i,j]
        end
      end
    end
  end
end

Sd2Init = Sd2
cbObj = sqrt.(Sd2)
O = sum(cbObj)
ObarInit = O/n

##### Subsection 2.5. Transferring grid points
stratcy = Array{Int64,2}(n,1)
# stratcy = Array{Int64,1}(n)
stratcy = strat0
TotTransf = 0
TotCycle = 0
for cycle = 1:maxcycle
  transfers = 0
  u = randperm(n)
  for t = u
    Delta = 0
    change = 0
    A = stratcy[t]
    ij = find(stratcy .== A)
    dA = sum(d2[t,ij])
    sumd2tinA = dA
    Sd2Amint = Sd2[A] - sumd2tinA
    cbObjA = sqrt.(abs.(Sd2Amint))
    for stratnr = 1:H
      Delta = 0
      sumd2plus = 0
      if stratnr != A
        B = stratnr
        ij = find(stratcy .== B)
        dB = sum(d2[t,ij])
        sumd2plus = dB
        cbObjB = sqrt.(abs.(Sd2[B] + sumd2plus))
        Delta = cbObjA + cbObjB - cbObj[A] -cbObj[B]
        if Delta < O*1e-10
          change = 1
          transfers = transfers + 1
          stratcy[t] = B            # update stratification
          Sd2[A] = Sd2[A] - sumd2tinA
          Sd2[B] = Sd2[B] + sumd2plus
          cbObj = sqrt.(abs.(Sd2))
          Obj = sum(cbObj)
          Delta = 0
        end                       # if Delta < Obj*1e-10
      end                         # if stratnr != A
      if change ==1 break end
    end                           # for strat=1:H
  end                             # for t=u

  # println("cycle ", cycle, "     transfers = ", transfers)
  TotTransf = TotTransf + transfers
  if transfers == 0 break end     # stopping rule
  TotCycle = cycle
end                               # for cycle=1:maxcycle

println("Total number of transfers = ", TotTransf)
println("Number of iteration cycles = ", TotCycle)

O = sum(cbObj)
Obar_sample = O/n
ObarFinal = Obar_sample

##### Subsection 2.6. Calculate size of sample-strata
Nh_sample = Array{Int64}(1,H)
Nh_sample = zeros(Int64,H,1)
for h = 1:H
  k = find(stratcy .== h)
  Nh_sample[h] = length(k)
end

##### SECTION 3. ALLOCATION TO THE SAMPLE-STRATA

# println("-------------------- Start allocating ----")

##### Subsection 3.1 Define the initial grid stratification of
##### length N, with zeros for the non-sample points,
##### and stratcy for the sample points.

stratification = zeros(Int64,N,1)
stratification[sample] = stratcy

sumsample = cbObj.^2   # sums of distances within sample-strata

##### Subsection 3.2 Partition the sample data according to the
##### sample-stratification.
Nh_sample = Array{Int64}(1,H)
deltas = Array{Float64,2}(N,H)
deltas = zeros(N,H)

for strat = 1:H
   k = find(stratification .== strat)
   Nh_sample[strat] = length(k)

##### Subsection 3.3 Extract a dataset for each variable from
##### sample-stratum "strat".
 xset = x[k]
 yset = y[k]
 z_predset = z_pred[k]
 s2set =s2[k]

##### Subsection 3.4 Calculate the sum of d2 for all gridpoints to
##### the points in sample-stratum "strat" using the datasets.
##### Calculate from the sum the increase of Obj (delta) if point
##### were allocated to "strat". Store deltas in N x H array
##### "deltas".
for i =1:N
 sumd2 = 0
 delta = 0
  for j = 1:Nh_sample[strat]
    d2ij = ((z_pred[i] - z_predset[j])^2)/R2 + (s2[i] + s2set[j]).*(1 - exp.(-3*(sqrt.((x[i] - xset[j])^2 + (y[i] - yset[j])^2)/range)))
    sumd2 = sumd2 + d2ij
  end    # j = 1:Nh_sample[strat]
 delta = sqrt(sumsample[strat] + sumd2) - cbObj[strat]
 deltas[i,strat] = delta
 end   # for i =1:N
end   # for strat = 1:H

##### Subsection 3.5 For each point t, find the stratum for which t
##### has the smallest delta, and allocate t to that stratum.
 for t = 1:N
  best = sortperm(deltas[t,1:H])
  stratification[t] = best[1]
 end    # for t = 1:N

#  println("---- End allocating ----")

##### SECTION 4 SAMPLE SIZES

##### Subsection 4.1 Size of the grid-strata.
Nh = Array{Int64}(1,H)
for strat = 1:H
  k = find(stratification .== strat)
  Nh[strat] = length(k)
end

##### Subsection 4.2 Contributions from grid-strata to O (as in
##### Subsect. 2.4). Distances between point-pairs in sample
##### are already calculated.
# println("- Calculating contribution from grid-strata to O -")

Sd2_grid = Array{Float64,1}(H)
Sd2_grid = zeros(Float64,H,1)
for i = 1:(N-1)
  for j = (i+1):N
    strati = stratification[i]
    stratj = stratification[j]
    if strati == stratj
      d2ij = ((z_pred[i] - z_pred[j])^2)/R2 + (s2[i] + s2[j]).*(1 - exp.(-3*(sqrt.((x[i] - x[j])^2 + (y[i] - y[j])^2)/range)))
      Sd2_grid[strati] = Sd2_grid[strati] + d2ij
    end
  end
end
Sd2_grid = Sd2_grid
cbObj_grid = sqrt.(Sd2_grid)
O = sum(cbObj_grid)
Obar_grid = O/N

##### Subsection 4.3 Total sample size before correction for
##### roundoff error.
n_pred = (CP*Area*Z_gamma*Obar_grid/(f*sqrt(2)))^(2/3)

##### Subsection 4.4 Neyman allocation.
sum_ahOh = 0
for h = 1:H
   ahOh = Nh[h]*cbObj_grid[h]
   sum_ahOh = sum_ahOh + ahOh
end

nh = Array{Float64}(1,H)
nhbest = Array{Float64}(1,H)
for h = 1:H
  nh[h] = n_pred *Nh[h] *cbObj_grid[h]/sum_ahOh
end
nh = round.(nh)

##### Subsection 4.5 Correction of total sample size to
##### avoid difference with sum of sizes in strata.
n_pred = sum(nh)
n_pred = convert(UInt64,n_pred)

##### Subsection 4.6 Check on smallest sample size in strata.
nh_low = indmin(nh)
nh_l = nh[nh_low]
nh_min = Array{Float64,1}
nh_min = nh_minim

##### Subsection 4.7 Update and output of intermediate results.
  Hbest = H
  stratbest = stratification
  nhbest = nh
  nbest = n_pred

  println("Intermediate results : " )
  # println("Obar : ", ObarFinal)
  println("Sample size :   ", nbest)
  println("Sample sizes in strata : ", nhbest)
  println("Smallest sample size allocated to a stratum :  ", nh_l)

if nh_l >= nh_min break end      # stop lowering H

end                                         # end optimization of H

##### SECTION 5. STRATIFIED RANDOM SAMPLING WITH FINAL
##### STRATIFICATION AND NEYMAN ALLOCATION

n_tot = sum(nhbest)
n_tot = round(n_tot)
n_tot = convert(UInt64,n_tot)
points = Array{UInt64}(1,n_tot)
xs = Array{Float64}(1,n_tot)
ys = Array{Float64}(1,n_tot)

for h = 1:Hbest
  k = find(stratbest .== h)      # numbers of points in stratum h
  stratsize = length(k)
  v = collect(1:stratsize)
  w = v[randperm(rng, stratsize)] # randomized indexes to points

  f=0
  for i = 1:nhbest[h]
    f=f+1                      # making Neyman allocations integer
  end

  k_rand = k[w]                # randomized points in stratum h
  points_h = k_rand[1:f]       # put the first f points in sample

  if h == 1                    # concatenate all H vectors of points
    points = points_h
  elseif h > 1
    points = vcat(points,points_h)
  end
end
xs = x[points]                 # get x coordinates of sample points
ys = y[points]                 # get y coordinates of sample points
strata = stratbest[points]   # get stratum numbers of sample points
sampnr = collect(1:n_tot)      # make sample numbers

strs = DataFrame()
strs[:SampleNr] = sampnr
strs[:StratNr] = strata
strs[:PointNr] = points
strs[:X] = xs
strs[:Y] = ys

stratXY = hcat(x,y,stratbest)

##### SECTION 6. OUTPUT OF FINAL RESULTS

writedlm("Stratification", stratXY)
CSV.write("Sample", strs)
println("  ")
println("FINAL RESULTS: ")
# println("Obar (O/N) : ", Obar_grid)
println("Number of strata : ", Hbest)
println("Total sample size : ", nbest)
println("Sample sizes in strata : ", nhbest)
println("Smallest sample size allocated to a stratum :  ", nh_l)
println("Size of grid-strata : ", Nh)

println("---- END OF FUNCTION OSPALL ---- ")

end   # function ospall
