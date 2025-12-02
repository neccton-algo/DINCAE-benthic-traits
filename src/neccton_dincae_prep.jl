# Prepare the files for DINCAE

using Pkg
Pkg.activate("neccton",shared=true)

using CSV
using DIVAnd
using DataFrames
using Glob
using Interpolations
using Missings
using NCDatasets
using PyPlot
using Random
using Statistics
using Dates
using DataStructures

                                        # Include("common.jl")

include("neccton_common.jl")
                # station_fname = joinpath(datadir,"Stations","stations.csv")
                # CWM_response_fname = joinpath(datadir,"Community Weighted Mean (sites x traits)","CWM_response.csv") N EXISTE PAS
                # env_matrix_fname = joinpath(datadir,"Environnement","matrice environnement.txt")
                # moddir = "/home/abarth/Data/NECCTON/PourSeverine"
                # figdir = expanduser("~/Figures/NECCTON")
                # split_fname = joinpath(basedir,"split.nc")


# Saving dataset function

function savedata(values,lon,lat,time,id,udates,varname,outname)
    len = length.(values);
    values = reduce(vcat,values);
    lon = reduce(vcat,lon);
    lat = reduce(vcat,lat);
    time = reduce(vcat,time);
    id = reduce(vcat,id);

    ds = NCDataset(outname,"c")

    defVar(ds,"size",len,("track",); attrib = OrderedDict(
        "sample_dimension" => "time"));
    defVar(ds,"dates",DateTime.(udates),("track",))

    defVar(ds,varname,values,("time",))
    defVar(ds,"lon",lon,("time",))
    defVar(ds,"lat",lat,("time",))
    defVar(ds,"id",id,("time",))
    defVar(ds,"dtime",time,("time",), attrib = OrderedDict(
        "long_name" => "time of measurement"))

    close(ds);
end



# Create df from station and CWM to get traits per location

df, resp, station = df_load(station_fname,CWM_response_fname);



# We check that we have all stations

env_matrix = CSV.read(env_matrix_fname,DataFrame)
env_matrix = rename(env_matrix, Symbol("Station ID") => :sta)
@assert Set(resp.sta) ⊆ Set(station.sta) # Error if the condition is not fullfilled
setdiff(resp.sta,env_matrix.sta)



# Load covar

fnames0 = sort(glob("Cl*_1d_*_*_btrc_T_*-*.nc",moddir))
fnames = sort(glob("Cl*_1d_*_*_grid_T_*-*.nc",moddir))


# Loading coordiantes of co-variables

lon1, lat1 = neccton_load_coord(fnames)


# Random definition of the trained part of the dataset

if !isfile(split_fname)
    Random.seed!(42)
    split_frac = 0.75
    p = randperm(size(df,1))
    nsplit = round(Int,split_frac * length(p))
    index_train = p[1:nsplit]
    index_val = p[nsplit+1:end]
    NCDataset(split_fname,"c") do ds
        defVar(ds,"index_train",index_train,("sample_train",))
        defVar(ds,"index_val",index_val,("sample_val",))
    end
end




# Separating training and validation

ds = NCDataset(split_fname)
index_train = ds["index_train"][:]
index_val = ds["index_val"][:]
close(ds)


# Bathymetry and mask loading

bathname = expanduser("~/Reconstruct_Points/Bathy/DivaData/Global/gebco_30sec_4.nc")
bathisglobal = true
mask,(pm,pn),(xi,yi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat)
hx, hy, h = DIVAnd.load_bath(bathname, bathisglobal, gridlon, gridlat)
#pcolormesh(hx,hy,h')

# Create a file for the mask

maskname = joinpath(basedir,"Results","mask.nc")
if !isfile(maskname)
    NCDataset(maskname,"c") do ds
        defVar(ds,"mask",Int8.(mask),("lon","lat"))
    end
end



# Longitude and latitude of train dataset

x = df.Longitude[index_train]
y = df.Latitude[index_train]


# Longitude and latitude of validation dataset

x_val = df.Longitude[index_val]
y_val = df.Latitude[index_val]


# Loading traits modalities

# Only Severine usefull traits
#n = "T15.M1"
#traits_modalities = ["T15.M1","T15.M2","T15.M3","T15.M4","T15.M5","T16.M1","T16.M2","T16.M3","T16.M4","T16.M5","T17.M1","T17.M2","T17.M3","T17.M4","T17.M5","T27.M1","T27.M2","T27.M3","T27.M4"]


traits_modalities = filter(s -> !isnothing(match(r"T.*\.M.*",s)),String.(propertynames(df)))
n = "T1.M1"


# For Analysis

for n in traits_modalities
    @show n
    @info "processing" n
    v = df[index_train,n]

    lon = [x]
    lat = [y]
    dtime = [fill(DateTime(1,1,1),size(x))]
    value = [v]
    id = [collect(1:length(v))]
    udates = [DateTime(1,1,1)]

    v = df[index_train,n]

    varname = n

    outname = joinpath(basedir,"DINCAE",varname * ".nc")
    @show outname
    mkpath(dirname(outname))

    savedata(value,lon,lat,dtime,id,udates,varname,outname)
end


# For validation

for n in traits_modalities
    @show n
    @info "processing" n
    v_val = df[index_val,n]

    lon_val = [x_val]
    lat_val = [y_val]
    dtime_val = [fill(DateTime(1,1,1),size(x_val))]
    value_val = [v_val]
    id_val = [collect(1:length(v_val))]
    udates = [DateTime(1,1,1)]

    v_val = df[index_val,n]

    varname = n

    outname = joinpath(basedir,"DINCAE",varname * "_val.nc")
    @show outname
    mkpath(dirname(outname))

    savedata(value_val,lon_val,lat_val,dtime_val,id_val,udates,varname,outname)
end



