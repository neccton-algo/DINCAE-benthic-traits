using DIVAnd
using CSV
using DataFrames
using PyPlot
using Statistics



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

function validate(fi)
    itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN)

    fi_val = itp.(df.Longitude[index_val],df.Latitude[index_val])


    v_val = df[index_val,n]

    @show sqrt(mean(filter(isfinite,(fi_val - v_val).^2)))
    @show sqrt(mean(filter(isfinite,(vm .- v_val).^2)))
    @show std(v)
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
pcolormesh(mask')


# Var Definition

# for(n qui tourne sur toutes les valeurs T)


n = "T1.M1" # Trait
v = df[index_train,n] 



#size(x)


# Parametrization

len = Float64(50000)
epsilon2 = 0.1


# DIVAndrun

# Compute anomaly
vm = mean(v)
v_anomaly = v.-vm

va,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),v_anomaly,len,epsilon2)


# Add mean to anomaly
va = va .+vm

pcolor(va')



# Validation

validate(va)


# save va + le resultat de validate
# end de toute les valeurs qui tournent sur t


