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

    close(ds)
end

include("neccton_common.jl")


station_fname = joinpath(datadir,"Stations","stations.csv")
CWM_response_fname = joinpath(datadir,"CWM_SxT","CWM_response.csv")
env_matrix_fname = joinpath(datadir,"Environnement","matrice environnement.txt")


#--


station = CSV.read(station_fname,DataFrame);

resp = CSV.read(CWM_response_fname,DataFrame,decimal=',');
resp = rename(resp, :Column1 => :sta)


env_matrix = CSV.read(env_matrix_fname,DataFrame)
env_matrix = rename(env_matrix, Symbol("Station ID") => :sta)

@assert Set(resp.sta) ⊆ Set(station.sta)

#@assert Set(resp.sta) ⊆ Set(env_matrix.sta)
#@assert Set(env_matrix.sta) ⊆ Set(station.sta)

setdiff(resp.sta,env_matrix.sta)

df = innerjoin(resp, station, on = :sta)


#plot(df.Longitude,df.Latitude,".")

figdir = expanduser("~/Figures/NECCTON")
#=
for n in names(df)
    d = df[:,n]
    if d[1] isa Number
        @show n
        clf();
        scatter(df.Longitude,df.Latitude,10,df[:,n])
        colorbar();
        title(n)
        savefig(joinpath(figdir,n * ".png"))
    end
end
=#

moddir = expanduser("~/Data/NECCTON/PourSeverine")
fnames0 = sort(glob("Cl*_1d_*_*_btrc_T_*-*.nc",moddir))

fnames = sort(glob("Cl1992_1d_*_*_grid_T_*-*.nc",moddir))


ds = NCDataset(fnames,"r")

lon = ds["nav_lon"][:,:]
lat = ds["nav_lat"][:,:]

lon = allowmissing(lon)
lon[lon .== -1] .= missing

lat = allowmissing(lat)
lat[lat .== -1] .= missing

lon1 = mapslices(s -> first(skipmissing(s)),lon,dims=2)[:,1]
lat1 = mapslices(s -> first(skipmissing(s)),lat,dims=1)[1,:]

#=
varname = "sbtemper"
sbtemper = nomissing(ds[varname][:,:,:],NaN)

mask = .!(isnan.(sbtemper[:,:,1]) .|| sbtemper[:,:,1] .== 0)

#pcolormesh(mask')

sbtemper[sbtemper .== 0] .= NaN

clf(); pcolormesh(lon1,lat1,sbtemper[:,:,1]')

mean_sbtemper = mean(sbtemper,dims=3)[:,:,1]
std_sbtemper = std(sbtemper,dims=3)[:,:,1]


ncovar = 2
field0 = zeros(length(lon1),length(lat1),ncovar)
field0[:,:,1] = mean_sbtemper
field0[:,:,2] = std_sbtemper

pcolormesh(lon1,lat1,std_sbtemper',cmap="jet")
=#

#pcolormesh(lon1,lat1,mean(sbtemper,dims=3)[:,:,1]'); colorbar()


#split between validation and training dataset

split_fname = joinpath(basedir,"split.nc")

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



bathname = expanduser("~/Data/DivaData/Global/gebco_30sec_4.nc")
bathisglobal = true

mask,(pm,pn),(xi,yi) = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat)
hx, hy, h = DIVAnd.load_bath(bathname, bathisglobal, gridlon, gridlat)

#pcolormesh(hx,hy,h')



x = df.Longitude[index_train]
y = df.Latitude[index_train]


traits_modalities = ["T1.M1","T1.M2","T1.M3","T1.M4","T1.M5"]

traits_modalities = filter(s -> !isnothing(match(r"T.*\.M.*",s)),String.(propertynames(df)))

n = "T1.M1"
for n in traits_modalities
    @info "processing" n
    v = df[index_train,n]

    positions = collect(zip(x,y))
    upositions = unique(positions)

    # lon = Vector{Float64}[]
    # lat = Vector{Float64}[]
    # dtime = Vector{DateTime}[]
    # value = Vector{Float64}[]
    # id = Vector{Int64}[]
    # udates = Vector{DateTime}(undef,length(upositions))

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

#=
for n in traits_modalities

len = 70e3
len = 80e3
epsilon2 = 0.5

vm = mean(v)
va = v .- vm;
fi,s = DIVAndrun(mask,(pm,pn),(xi,yi),(x,y),va,len,epsilon2);
fi .+= vm

cl = quantile(df[:,n],(0.1,0.9))

function plmap(;orientation = "horizontal")
    clim(cl)
    colorbar(orientation=orientation);
    xlim(gridlon[[1,end]]); ylim(gridlat[[1,end]])
    plotmap();
    set_aspect_ratio()
end

#=
clf();
subplot(2,2,1); pcolor(xi,yi,fi); plmap()
subplot(2,2,2); scatter(x,y,10,v); plmap();
subplot(2,2,3); scatter(df.Longitude[index_val],df.Latitude[index_val],10,df[index_val,n]);
plmap()
=#


function validate(fi)
    itp = extrapolate(interpolate((gridlon,gridlat),fi,Gridded(Linear())),NaN)

    fi_val = itp.(df.Longitude[index_val],df.Latitude[index_val])


    v_val = df[index_val,n]

    @show sqrt(mean(filter(isfinite,(fi_val - v_val).^2)))
    @show sqrt(mean(filter(isfinite,(vm .- v_val).^2)))
    @show std(v)
end

validate(fi)





=#
