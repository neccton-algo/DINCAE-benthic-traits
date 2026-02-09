using DIVAnd
using CSV
using DataFrames
using PyPlot
using Statistics
using Glob
using JSON3
using OceanPlot

# Common.jl

include("neccton_common.jl")

# Load Dataset

df, resp, station = df_load(station_fname,CWM_response_fname);

# Split dataset 

if !isfile(split_fname)
    Random.seed!(42)
    split_frac = 0.75
    p = randperm(size(df,1))
	println(size(p[:]))
    nsplit = round(Int,split_frac * length(p))
    index_train = p[1:nsplit]
    index_val = p[nsplit+1:end]
    NCDataset(split_fname,"c") do ds
        defVar(ds,"index_train",index_train,("sample_train",))
        defVar(ds,"index_val",index_val,("sample_val",))
    end
end


# Validate

function validateDIVA(fi,n) # fi = DIVAnd result, n = varname

    itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN) # Creation of the interpolator of the result fi

    fi_val = itp.(df.Longitude[index_val],df.Latitude[index_val]) # On the validation data


    v_val = df[index_val,n] # Taking values of the validation dataset on the considered stations

    RMS1 = sqrt(mean(filter(isfinite,(fi_val - v_val).^2))) 
        
    return RMS1
    
end


function plmap(cl;orientation = "horizontal")
    clim(cl)
    colorbar(orientation=orientation);
    xlim(gridlon[[1,end]]); ylim(gridlat[[1,end]])
    OceanPlot.plotmap();
    OceanPlot.set_aspect_ratio()
end



# Separating training and validation following Split

ds = NCDataset(split_fname)
index_train = ds["index_train"][:]
index_val = ds["index_val"][:]
close(ds)


# Longitude and latitude

x = df.Longitude[index_train]
y = df.Latitude[index_train]


# Bathymetry and mask loading

bathname = expanduser("~/Reconstruct_Points/Bathy/DivaData/Global/combined_emodnet_bathymetry.nc")
bathisglobal = true
mask,(pm,pn),(xi,yi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat)
hx, hy, h = DIVAnd.load_bath(bathname, bathisglobal, gridlon, gridlat)
#pcolormesh(hx,hy,h')
#pcolormesh(mask')



# Parametrization

#=
len = Float64(12^4)
epsilon2 = 0.007
=#

len = 80000.0
epsilon2 = 1.

# Var Definition

#varname = "T2.M1" # Trait
#v = df[index_train,varname] 



varname = "T1.M1"
varnames = replace.(
    basename.(
        filter(f->!endswith(f,"_val.nc"), glob("T*M*.nc",joinpath(basedir,"DINCAE")))
        ),".nc" => "")


for varname in varnames

# Train set load
v = df[index_train,varname] # 158 valeurs de trainset


# DIVAndrun
cl = quantile(v[:],(0.001,0.90))

# With anomaly
vm = mean(v)
v_anomaly = v.-vm

va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v_anomaly,len,epsilon2)

va = va .+vm # Add mean to anomaly

fig1 = pcolor(xi,yi,va);plmap(cl);title("Baseline DIVAnd_anmoaly for $varname")


RMS_DIVAnd_anomaly = validateDIVA(va,varname)

    
# Without anomalies
vav,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,len,epsilon2)


pcolor(xi,yi,vav);plmap(cl);title("Basline DIVAnd for $varname")
savefig(joinpath(figdir, "DIVAnd-result-$varname.png"))

RMS_DIVAnd = validateDIVA(vav,varname)


    
# Save results
open("sortieDIVAnd.txt", "a") do f


    write(f,varname * "\n")
    
    JSON3.pretty(f,"RMS_anomaly = $(RMS_DIVAnd_anomaly)"; allow_inf=true)
    write(f, "\n")
    
    JSON3.pretty(f,"RMS = $(RMS_DIVAnd)"; allow_inf=true)

    write(f, "\n")

end
    


end


# save va + le resultat de validate
# end de toute les valeurs qui tournent sur t
