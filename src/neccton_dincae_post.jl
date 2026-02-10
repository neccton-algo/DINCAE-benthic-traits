# Post-processing for the DINCAE results

using CSV
using DataFrames
using OceanPlot
using PyPlot, Statistics
using Interpolations
using Glob
using JSON3
using Dates
using FileWatching


# Validation function

function validate(n,fi,fi_err,epochs,probability_skip_for_training,learning_rate_decay_epoch,regularization_L2_beta,upsampling_method,learning_rate)


    fi_itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN) # Create fi interpolator
    fi_val = fi_itp.(df.Longitude[index_val],df.Latitude[index_val]) # Interpolate on validation index

    fi_train_itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN) # Create fi train interpolator
    fi_train = fi_train_itp.(df.Longitude[index_train], df.Latitude[index_train]) # Inteprolate fi on training points
    
    fi_err_itp = extrapolate(interpolate((gridlon,gridlat),fi, Gridded(Linear())),NaN) # Create fi_error interpolator
    fi_err_val = fi_err_itp.(df.Longitude[index_val],df.Latitude[index_val]) # Interpolate error on validation


    v_train = df[index_train,n] # Values of training
    v_val = df[index_val,n] # Values of validation
    sel = isfinite.(fi_val) .& isfinite.(v_val)
    sel_train = isfinite.(fi_train)
    
    sigmas = (;  (Symbol("sigma$i") => mean(abs.(v_val[sel]-fi_val[sel]) .<  i*fi_err_val[sel]) for i = 1:3)...)

    RMS_train = sqrt(mean((fi_train[sel_train] - v_train[sel_train]).^2))
    RMS = sqrt(mean((fi_val[sel] - v_val[sel]).^2))
    correlation = cor(fi_val[sel],v_val[sel])
    bias = mean(fi_val[sel] - v_val[sel])

    std_obs = std(v_val[sel])
    std_fi = std(fi_val[sel])

    epochs = epochs
    
    #@show RMS
    #    @show sqrt(mean(filter(isfinite,(vm .- v_val).^2)))


    return (; epochs,probability_skip_for_training,learning_rate_decay_epoch,regularization_L2_beta,upsampling_method,learning_rate, RMS, RMS_train, correlation, bias, std_obs, std_fi, sigmas...)
end


# Plot function

function plmap(cl;orientation = "horizontal")
    clim(cl)
    colorbar(orientation=orientation);
    xlim(gridlon[[1,end]]); ylim(gridlat[[1,end]])
    OceanPlot.plotmap();
    OceanPlot.set_aspect_ratio()
end



# Include prep to retrieve df variable

#include("neccton_dincae_prep.jl")


#figdir = joinpath(basedir,"Figures")
#mkpath(figdir)


# For any trait

varname = "T1.M1"
varnames = replace.(
    basename.(
        filter(f->!endswith(f,"_val.nc"), glob("T*M*.nc",joinpath(basedir,"DINCAE")))
        ),".nc" => "")

# Subset for first RMS analysis 1=>18 
#varnames = varnames[1:51]
#varnames = varnames[51:end]

# For severine exdidataset

#varname = "T15.M1"
#varnames = ["T15.M1","T15.M2","T15.M3","T15.M4","T15.M5","T16.M1","T16.M2","T16.M3","T16.M4","T16.M5","T17.M1","T17.M2","T17.M3","T17.M4","T17.M5","T27.M1","T27.M2","T27.M3","T27.M4"]

#time_run = "DINCAE_2025-02-25T143624"


# Grid

xi = gridlon .+ 0 * gridlat'
yi = 0 * gridlon .+ gridlat'


# Corrections suite au message d'erreur avec les ARGS

df, resp, station = df_load(station_fname,CWM_response_fname);
ds = NCDataset(split_fname)
index_train = ds["index_train"][:]
index_val = ds["index_val"][:]
close(ds)

x = df.Longitude[index_train]
y = df.Latitude[index_train]

xval = df.Longitude[index_val]
yval = df.Latitude[index_val]




######## For parrallel #########

# Index and Args

#index = Int[]
#index = parse(Int,ARGS[1])
# replace varname by the index of the sbatch
# Each sbatch iteration change varname
#varname = varnames[index]


# Lock file for parralelisation
#FileWatching.Pidfile.mkpidlock(expanduser("~/Reconstruct_Points/DINCAE-benthic-traits/src/DINCAE.pid")) do

#################################

varname = "T1.M1"
n = varname


dsmodel = NCDataset(joinpath(auxdir,"resized_clim_2008_2018.nc"))
avg_ox = dsmodel["avg_oxygen"][:,:]
close(dsmodel)

mask = avg_ox .>= 5


# Seeking for the varnames in the dataset
open(joinpath(outdir,"sortieDINCAE.txt"), "a") do f

#for varname in ["T1.M1"]

    local ds, fnames_rec, v, fi, fi_err, n, cl

    @show varname
    
    v = df[index_train,varname]
    vval = df[index_val,varname]
    
    fnames_rec = [joinpath(outdir,"data-avg-$varname.nc")]
    ds = NCDataset(fnames_rec[1])

    # Mean
    dsanom = NCDataset(
    joinpath(basedir, "DINCAE-anom",varname*"-anom.nc"),"a")

    fi = nomissing(ds[varname][:,:,1])
    fi = fi .+ dsanom.attrib["mean"]
    
    fi_err = nomissing(ds[varname * "_error"][:,:,1])
    fi_err = fi_err .+ dsanom.attrib["mean"]
    

    close(dsanom)

    n = varname
    
    # validate function on the current variable
    summary = validate(n,fi,fi_err,epochs,probability_skip_for_training,learning_rate_decay_epoch,regularization_L2_beta,upsampling_method,learning_rate)
    @show summary

    
# Plot + validation statistics computation

      
    n = varname
    # validate function on the current variable
    summary = validate(n,fi,fi_err,epochs,probability_skip_for_training,learning_rate_decay_epoch,regularization_L2_beta,upsampling_method,learning_rate)
    @show summary

    # apply mask only for the plot
    fi[.!mask] .= NaN
    fi_err[.!mask] .= NaN
    
# Plot + validation statistics computation
    cl = quantile(df[:,n],(0.1,0.9))

    fi[.!mask] .= NaN
    fi_err[.!mask] .= NaN
    
    clf();
    subplot(1,4,1); scatter(x,y,10,v); plmap(cl); title("Observation ($n)")
    subplot(1,4,1); scatter(xval,yval,10,vval); plmap(cl); title("Validation ($n)")
    subplot(1,4,2); pcolor(xi,yi,fi); plmap(cl); title("Analysis")
    subplot(1,4,3); pcolor(xi,yi,fi_err); plmap((0, 0.16)); title("std. err.")
    savefig(joinpath(figdir,"analysis-$n.png"))

    # Writing results in f "sortieDINCAE"
    
    write(f,varname * "\n")   
    JSON3.pretty(f, summary; allow_inf=true)
    write(f, "\n")
    
    end
 #    end # stop writing

#end # close lock






