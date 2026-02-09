# Execute DINCAE

using Dates
using DINCAE
using LinearAlgebra
using NCDatasets
using Printf
using Random
using DIVAnd
using CUDA
#using DINCAE_altimetry
using Glob
using JSON3
using DataStructures
using IntervalSets


#include("neccton_dincae_prep.jl")
include("neccton_common.jl")

# Set the precision : = Float32 for fast run
const T = Float32
const F = Float32

# Resolution coming from common.jl

Δlon = dlon
Δlat = dlat

lonr = gridlon
latr = gridlat

@show size(lonr)




# Parametrization of DINCAE

# Modified ones

#epochs = 1000
epochs = 1000
#epochs = 1500


#probability_skip_for_training = 0.6 AllCovar
probability_skip_for_training = 0.1

learning_rate_decay_epoch = 500

regularization_L2_beta = 0

#upsampling_method = :nearest
upsampling_method = :bilinear

learning_rate = 0.004906558111668519
# 0.000490 AllCovar

# Unmodified ones

batch_size = 1;

clip_grad = 5.0


#save_epochs = 60:10:epochs
save_epochs = [epochs]

#enc_nfilter_internal = [25,50,75]
enc_nfilter_internal = round.(Int,32 * 2 .^ (0:5))


#jitter_std_pos = 0.5 .* (0.17145703272237467f0,0.17145703272237467f0)
jitter_std_pos = (0f0,0f0)

ntime_win = 1



skipconnections = 3:length(enc_nfilter_internal)
#skipconnections = 3:(length(enc_nfilter_internal)+1)

laplacian_penalty = 0.8f0
laplacian_penalty = 0.9f0


loss_weights_refine = (1,)
#loss_weights_refine = (0.1,0.2,0.7)

seed = 12345


# Parametrization of parameters that have to be random

#learning_rate = Float32(10 ^ (rand(-4..(-3))))
#laplacian_penalty = Float32(10 ^ (rand(-0..1)))
#enc_nfilter_internal = round.(Int,32 * 2 .^ (0:rand(2:5)))

#regularization_L2_beta = Float32(10 ^ (rand(-6..(-1))))
#epochs = rand(300:1500)
#upsampling_method = rand([:nearest,:bilinear])
#save_epochs = 200:rand(20:60):epochs
save_epochs = [epochs]

#epochs = 300



# Temporal scale definition

timestamp = Dates.format(Dates.now(),"yyyy-mm-ddTHHMMSS")
outdir = joinpath(basedir,"Results","DINCAE-$(timestamp)")
mkpath(outdir)
paramsname = joinpath(outdir,"params.json")


# Defining the grid based on the resolution

grid = (lonr,latr)



# GPU availability check

if CUDA.functional()
    Atype = CuArray{F}
else
    @warn "No supported GPU found. We will use the CPU which is very slow. Please check https://developer.nvidia.com/cuda-gpus"
    Atype = Array{F}
end

@show Atype

#cp(@__FILE__,joinpath(outdir,basename(@__FILE__)),force=true)


# load covariables

covars_fname = [
    # (filename = "mean_oxygenbot.nc", varname = "mean_oxygenbot", errvarname = nothing),
    # (filename = "std_oxygenbot.nc",  varname = "std_oxygenbot",  errvarname = nothing),
    # (filename = "mean_DOX.nc",      varname = "mean_dox",      errvarname = nothing),
    # (filename = "std_DOX.nc",       varname = "std_dox",       errvarname = nothing),
    # (filename = "low_DOX.nc",       varname = "low_dox",       errvarname = nothing),
    (filename = "resized_clim_2008_2018.nc",       varname = "avg_oxygen",       errvarname = nothing),
    (filename = "resized_clim_2008_2018.nc",       varname = "avg_fCSED",       errvarname = nothing),
    (filename = "resized_clim_2008_2018.nc",       varname = "avg_sCSED",       errvarname = nothing),
    (filename = "resized_clim_2008_2018.nc",       varname = "avg_POC",       errvarname = nothing),
    (filename = "resized_clim_2008_2018.nc",       varname = "avg_Botflux",       errvarname = nothing),
    (filename = "resized_clim_2008_2018.nc",       varname = "Bath",       errvarname = nothing),

    (filename = "resized_sediments2.nc",       varname = "sediment_type1",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type2",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type3",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type4",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type5",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type6",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type7",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type8",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type9",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type10",       errvarname = nothing),
    (filename = "resized_sediments2.nc",       varname = "sediment_type11",       errvarname = nothing),
]

# Put everything in auxdata_files

auxdata_files = map(entry -> (;entry...,filename = joinpath(auxdir,entry.filename)),covars_fname)



# Variables definition

varname = "T1.M1"
varnames = replace.(
    basename.(
        filter(f->!endswith(f,"_val.nc"),glob("T*M*.nc",joinpath(basedir,"DINCAE")))
        ),".nc" => "")




# Choice of the variables

#varnames = varnames[1:51]
#varnames = varnames[52:end]

# For ARGS use: do a loop on sbatch which will take varnames[index] for variable
#for i in $(seq 118),   ; do sbatch neccton_dincae_run.jl $i 
# index = 1
# Index being a number and not a string
#index = parse(Int,ARGS[1])
# replace varname by the index of the sbatch
index = 1
varname = varnames[index]

mkpath(outdir)






# DINCAE.reconstruct_points for every variable

#for varname in varnames
    filename = joinpath(basedir,"DINCAE",varname * ".nc")

    fnames_rec = [joinpath(outdir,"data-avg-$varname.nc")] # fichier de sauvegarde au sein de rp.


# Parameter saved in a directory

open(paramsname,"w") do f
    JSON3.pretty(f,OrderedDict(
        "laplacian_penalty" => laplacian_penalty,
        "epochs" => epochs,
        "batch_size" => batch_size,
        "enc_nfilter_internal" => enc_nfilter_internal,
        "clip_grad" => clip_grad,
        "regularization_L2_beta" => regularization_L2_beta,
        "ntime_win" => ntime_win,
        "upsampling_method" => upsampling_method,
        "loss_weights_refine" => loss_weights_refine,
        "save_epochs" => save_epochs,
        "skipconnections" => skipconnections,
        "dlon" => dlon,
        "dlat" => dlat,
        "jitter_std_pos" => jitter_std_pos,
        "probability_skip_for_training" => probability_skip_for_training,
        "learning_rate" => learning_rate,
        "learning_rate_decay_epoch" => learning_rate_decay_epoch,
        "covars_fname" => covars_fname,
        "seed" => seed,
    ))
end


DINCAE.reconstruct_points(
    T,Atype,filename,varname,grid,fnames_rec;
    learning_rate = learning_rate,
    learning_rate_decay_epoch = learning_rate_decay_epoch,
    epochs = epochs,
    batch_size = batch_size,
    enc_nfilter_internal = enc_nfilter_internal,
    skipconnections = skipconnections,
    clip_grad = clip_grad,
    save_epochs = save_epochs,
    upsampling_method = upsampling_method,
    jitter_std_pos = jitter_std_pos,
    probability_skip_for_training = probability_skip_for_training,
    ntime_win = ntime_win,
    laplacian_penalty = laplacian_penalty,
    auxdata_files = auxdata_files,
    loss_weights_refine = loss_weights_refine,
)
#end

 
# Validation through dincae_post


include("neccton_dincae_post.jl")




