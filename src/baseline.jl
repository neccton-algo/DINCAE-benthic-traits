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

function validateDIVA(fi,n)
    itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN) # Interpolation of the result

    fi_val = itp.(df.Longitude[index_val],df.Latitude[index_val]) # On the validation data


    v_val = df[index_val,n] # Taking the variable name

    RMS1 = sqrt(mean(filter(isfinite,(fi_val - v_val).^2))) 
    #RMS2 = sqrt(mean(filter(isfinite,(vm .- v_val).^2)))
    

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

len = Float64(10^4)
epsilon2 = 0.01



# Var Definition

#varname = "T2.M1" # Trait
#v = df[index_train,varname] 



varname = "T1.M1"
varnames = replace.(
    basename.(
        filter(f->!endswith(f,"_val.nc"), glob("T*M*.nc",joinpath(basedir,"DINCAE")))
        ),".nc" => "")


for varname in varnames

    v = df[index_train,varname] # 158 valeurs de trainset


# DIVAndrun

cl = quantile(v[:],(0.01,0.99))

# Compute anomaly
vm = mean(v)
v_anomaly = v.-vm

va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v_anomaly,len,epsilon2)

# Add mean to anomaly
va = va .+vm

pcolor(va')
colorbar()

# Without anomalies
vav,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v,len,epsilon2)

   fig =  pcolor(xi,yi,vav);plmap(cl);title("Basline DIVAnd for $varname")
    savefig(fig)
    

RMS_DIVAnd = validateDIVA(va,varname)
    
# Validation
open("sortieDIVAnd.txt", "a") do f


    write(f,varname * "\n")
    
    JSON3.pretty(f,"RMS = $(RMS_DIVAnd)"; allow_inf=true)
    write(f, "\n")

end
    


end


# save va + le resultat de validate
# end de toute les valeurs qui tournent sur t
