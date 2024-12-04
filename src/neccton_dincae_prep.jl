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



                                        # Reading and dataset loading with common names


df, resp, station = df_load(station_fname,CWM_response_fname);


env_matrix = CSV.read(env_matrix_fname,DataFrame)
env_matrix = rename(env_matrix, Symbol("Station ID") => :sta)
@assert Set(resp.sta) ⊆ Set(station.sta) #verifie que la condition est vraie sinon erreure
# on regarde que toutes las stations dans response sont cpmprise dans sta


setdiff(resp.sta,env_matrix.sta)



                                    # Load covar

fnames0 = sort(glob("Cl*_1d_*_*_btrc_T_*-*.nc",moddir))
fnames = sort(glob("Cl*_1d_*_*_grid_T_*-*.nc",moddir))


                                    # Function defined in common


lon1, lat1 = neccton_load_coord(fnames)


                                  # Random definition of the trained part of the dataset
                                  # Same lines in post
                                


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


ds = NCDataset(split_fname)
index_train = ds["index_train"][:]
index_val = ds["index_val"][:]
close(ds)

                                                # Bathymetry and mask loading

bathname = expanduser("~/Data/DivaData/Global/gebco_30sec_4.nc")
bathisglobal = true
mask,(pm,pn),(xi,yi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat)
hx, hy, h = DIVAnd.load_bath(bathname, bathisglobal, gridlon, gridlat)
#pcolormesh(hx,hy,h')


x = df.Longitude[index_train]
y = df.Latitude[index_train]

                                     # Definition of the modalitites that are trained + save


#traits_modalities = ["T1.M1","T1.M2","T1.M3","T1.M4","T1.M5"]

traits_modalities = filter(s -> !isnothing(match(r"T.*\.M.*",s)),String.(propertynames(df)))

n = "T1.M1"
for n in traits_modalities
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



