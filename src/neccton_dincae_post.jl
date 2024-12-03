                                                # Post-processing for the DINCAE results


using CSV
using DataFrames
using PyPlot, Statistics
using OceanPlot
using Interpolations
using Glob
using JSON3


                                                       # Validation function

function validate(n,fi,fi_err)
    fi_itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN)
    fi_val = fi_itp.(df.Longitude[index_val],df.Latitude[index_val])

    fi_err_itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN)
    fi_err_val = fi_err_itp.(df.Longitude[index_val],df.Latitude[index_val])

    v_val = df[index_val,n]
    sel = isfinite.(fi_val) .& isfinite.(v_val)

    sigmas = (;  (Symbol("sigma$i") => mean(abs.(v_val[sel]-fi_val[sel]) .<  i*fi_err_val[sel]) for i = 1:3)...)

    RMS = sqrt(mean((fi_val[sel] - v_val[sel]).^2))
    correlation = cor(fi_val[sel],v_val[sel])
    bias = mean(fi_val[sel] - v_val[sel])

    std_obs = std(v_val[sel])
    std_fi = std(fi_val[sel])
#    @show RMS
    #    @show sqrt(mean(filter(isfinite,(vm .- v_val).^2)))
    return (; RMS, correlation, bias, std_obs, std_fi, sigmas...)
end


                                                           # Plot function

function plmap(cl;orientation = "horizontal")
    clim(cl)
    colorbar(orientation=orientation);
    xlim(gridlon[[1,end]]); ylim(gridlat[[1,end]])
    OceanPlot.plotmap();
    OceanPlot.set_aspect_ratio()
end


                                                        # Loading of split_fname
                                          # The file that differentiate the training and the dataset 


include("neccton_common.jl")

df = df_load(station_fname,CWM_response_fname)

##### Repetition dincae_prep

#ds = NCDataset(split_fname)
#index_train = ds["index_train"][:]
#index_val = ds["index_val"][:]
#close(ds)


#x = df.Longitude[index_train]
#y = df.Latitude[index_train]

######


# outdir = joinpath(basedir,"DINCAE-temp")
# outdir = joinpath(basedir,"DINCAE-temp3-rerun")
# outdir = joinpath(basedir,"DINCAE-temp4")
# outdir = joinpath(basedir,"DINCAE-temp-all")
#outdir = joinpath(basedir,"DINCAE-temp-all-$(1/dlon)-3")
#outdir = "/home/abarth/Data/NECCTON/DINCAE//DINCAE-2024-02-23T161128/"

figdir = joinpath(outdir,"Fig")
mkpath(figdir)

varname = "T1.M1"
varnames = replace.(basename.(glob("T*M*.nc",joinpath(basedir,"DINCAE"))),".nc" => "")

xi = gridlon .+ 0 * gridlat'
yi = 0 * gridlon .+ gridlat'




                                              # Seeking for the varnames in the dataset

for varname in  varnames
    local ds, fnames_rec, v, fi, fi_err, n, cl

    v = df[index_train,varname]

    fnames_rec = [joinpath(outdir,"data-avg-$varname.nc")]

    ds = NCDataset(fnames_rec[1])
    fi = nomissing(ds[varname][:,:,1])
    fi_err = nomissing(ds[varname * "_error"][:,:,1])

  
    
    
    #=
    clf();
    subplot(2,2,1); scatter(x,y,10,v); plmap(); title("Observation, training ($n)")
    subplot(2,2,2); scatter(df.Longitude[index_val],df.Latitude[index_val],10,df[index_val,n]); title("Observation, validation ($n)")
    plmap()
    subplot(2,2,3); pcolor(xi,yi,fi); plmap()
    title("Analysis (using training data)")
    savefig(joinpath(figdir,"analysis-$n.png"))
    =#

                                               # Plot + validation statistics computation
    
    
    
    n = varname
    cl = quantile(df[:,n],(0.1,0.9))

    clf();
    subplot(1,3,1); scatter(x,y,10,v); plmap(cl); title("Observation ($n)")
    subplot(1,3,2); pcolor(xi,yi,fi); plmap(cl); title("Analysis")
    subplot(1,3,3); pcolor(xi,yi,fi_err); plmap((0, 0.16)); title("std. err.")
    savefig(joinpath(figdir,"analysis-$n.png"))

    summary = validate(n,fi,fi_err)
    @show summary
    statname = replace(fnames_rec[1],".nc" => ".json")
    open(statname,"w") do f
        JSON3.pretty(f,summary; allow_inf = true)
    end
end


