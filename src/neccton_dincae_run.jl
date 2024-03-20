# Execute DINCAE


using Pkg; Pkg.activate("/home/abarth/src/DINCAE-benthic-traits")

using Dates
using DINCAE
using LinearAlgebra
using NCDatasets
using Printf
using Random
using DIVAnd
using CUDA
using Glob
using JSON3
using DataStructures
using IntervalSets

include("neccton_common.jl")


T = Float32

bathname = expanduser("~/Data/DivaData/Global/gebco_30sec_4.nc")
bathisglobal = true

mask,(pm,pn),(xi,yi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat)

maskname = joinpath(basedir,"DINCAE","mask.nc")
if !isfile(maskname)
    NCDataset(maskname,"c") do ds
        defVar(ds,"mask",Int8.(mask),("lon","lat"))
    end
end

batch_size = 1;
Δlon = dlon
Δlat = dlat
lonr = gridlon
latr = gridlat


clip_grad = 5.0
epochs = 1000
#epochs = 500
save_epochs = 60:10:epochs
save_epochs = [epochs]
#enc_nfilter_internal = [25,50,75]
enc_nfilter_internal = [25,50,75,75,75,75]
enc_nfilter_internal = round.(Int,32 * 2 .^ (0:4))
upsampling_method = :nearest
#upsampling_method = :bilinear
probability_skip_for_training = 0.4
#probability_skip_for_training = 1.
jitter_std_pos = (0.17145703272237467f0,0.17145703272237467f0)
#jitter_std_pos = 0.5 .* (0.17145703272237467f0,0.17145703272237467f0)
jitter_std_pos = (0f0,0f0)
ntime_win = 1
learning_rate = 0.0004906558111668519
skipconnections = 3:length(enc_nfilter_internal)
skipconnections = 3:(length(enc_nfilter_internal)+1)
learning_rate_decay_epoch = 50
laplacian_penalty = 1f-3
laplacian_penalty = 1f-2
laplacian_penalty = 1f-1
laplacian_penalty = 2f-1
laplacian_penalty = 1f0
laplacian_penalty = 0.8f0
#laplacian_penalty = 0
regularization_L2_beta = 0
loss_weights_refine = (1,)
#loss_weights_refine = (0.1,0.2,0.7)
seed = 12345


# random
learning_rate = Float32(10 ^ (rand(-4..(-3))))
laplacian_penalty = Float32(10 ^ (rand(-0..1)))
enc_nfilter_internal = round.(Int,32 * 2 .^ (0:rand(2:5)))
#regularization_L2_beta = Float32(10 ^ (rand(-6..(-1))))
epochs = rand(300:1500)
upsampling_method = rand([:nearest,:bilinear])
#save_epochs = 200:rand(20:60):epochs
save_epochs = [epochs]

#epochs = 300
#------

timestamp = Dates.format(Dates.now(),"yyyy-mm-ddTHHMMSS")

outdir = joinpath(basedir,"DINCAE","DINCAE-$(timestamp)")

mkpath(outdir)

paramsname = joinpath(outdir,"params.json")

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
        "batch_size" => batch_size,
        "seed" => seed,
    ))
end

const F = Float32

grid = (lonr,latr)

if CUDA.functional()
    Atype = CuArray{F}
else
    @warn "No supported GPU found. We will use the CPU which is very slow. Please check https://developer.nvidia.com/cuda-gpus"
    Atype = Array{F}
end


cd(joinpath(dirname(pathof(DINCAE)),"..")) do
    if isdir(".git")
        write("$outdir/DINCAE.commit", read(`git rev-parse HEAD`))
        write("$outdir/DINCAE.diff", read(`git diff`))
    end
end;

if isfile(@__FILE__)
    cp(@__FILE__,joinpath(outdir,basename(@__FILE__)),force=true)
end

# function cvrms(fname_rec)
#     varname = "sla"

#     filename_dev = expanduser("~/tmp/Altimetry/all-sla.dev.nc")
#     fnamesummary_dev = replace(fname_rec,".nc" => ".dev.json")

#     filename_test = expanduser("~/tmp/Altimetry/all-sla.test.nc")
#     fnamesummary_test = replace(fname_rec,".nc" => ".test.json")

#     summary, = DINCAE_altimetry.errstat(filename_dev,fname_rec,varname,maskname; fnamesummary = fnamesummary_dev)
#     return summary["cvrms"]
# end

# load covariables
covars_fname = [
    (filename = "mean_sbtemper.nc", varname = "mean_sbtemper", errvarname = nothing),
    (filename = "std_sbtemper.nc",  varname = "std_sbtemper",  errvarname = nothing),
    (filename = "mean_DOX.nc",      varname = "mean_dox",      errvarname = nothing),
    (filename = "std_DOX.nc",       varname = "std_dox",       errvarname = nothing),
    (filename = "low_DOX.nc",       varname = "low_dox",       errvarname = nothing),
    (filename = "mean_PAR.nc",      varname = "mean_par",      errvarname = nothing),
    (filename = "std_PAR.nc",       varname = "std_par",       errvarname = nothing),
]

auxdata_files = map(entry -> (;entry...,filename = joinpath(auxdir,entry.filename)),covars_fname)


varname = "T1.M1"
varnames = replace.(basename.(glob("T*M*.nc",joinpath(basedir,"DINCAE"))),".nc" => "")

mkpath(outdir)
#for varname in varnames[1:1]
for varname in varnames
    filename = joinpath(basedir,"DINCAE",varname * ".nc")

    fnames_rec = [joinpath(outdir,"data-avg-$varname.nc")]


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
end

#@show cvrms(fnames_rec[1])
#end

include("neccton_dincae_post.jl")
