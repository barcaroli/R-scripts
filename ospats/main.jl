
# To start the calculations: run this script, after filling in the input
# filename and the process parameters below.
# Date: 2018-03-30
# Version: 1.01
# Author: Jaap de Gruijter
cd("C:\\Users\\UTENTE\\Google Drive\\Spatial Sampling\\meuse\\ospats")
using DataFrames
using CSV
include("readdata.jl")
include("ospats.jl")
include("ospall.jl")
filename = "example.txt"     ########### INSERT here name of data file:

println("=======================================================")
println("File name :  ",filename)

################################## INSERT here process parameters:
H_min = 3       # minimum number of strata.
H_max = 3     # maximum number of strata.
nh_minim = 2    # minimum sample size allowed in the strata.
CP = 1         # carbon offset price (A$ per Mg).
f = 1         # cost of obtaining data per grid point (A$).
Area = 648    # surface area of the farm (ha).
# Z_gamma = 1.645 # 90% quantile of the standard normal distribution.
Z_gamma = 1.96 # 95% quantile of the standard normal distribution.
# Z_gamma = 2.326 # 98% quantile of the standard normal distribution.
# Z_gamma = 2.576 # 99% quantile of the standard normal distribution.

R2 = 1       # squared correlation coefficient of predictions.
range = 520     # auto-correlation range of the prediction error,
                # in the same length unit as x and y.
maxcycle = 150  # Maximal number of re-allocation cycles.
in = 1           # Sampling interval for coarse-gridding
                # (in = 1 implies no coarse-gridding).
seed = 1234     # seed for the rand() and randperm() functions.

rng = srand(seed)
readdata()
if in == 1
  println("Calling function ospats")
  ospats()
elseif in > 1
  println("Calling function ospall")
  ospall()
end     # if in == 1
