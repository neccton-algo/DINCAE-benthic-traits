using NCDatasets
using Statistics
using Missings
using DIVAnd
using Glob

include("neccton_common.jl")

mkpath(auxdir)

fnames0 = sort(glob("Cl*_1d_*_*_btrc_T_*-*.nc",moddir))
fnames = sort(glob("Cl1992_1d_*_*_grid_T_*-*.nc",moddir))

fnames_Btrc = sort(glob("Cl1992_1d_*_*_Btrc_T_*-*.nc",moddir))

ds = NCDataset(fnames,"r")

lon = ds["nav_lon"][:,:]
lat = ds["nav_lat"][:,:]

lon = allowmissing(lon)
lon[lon .== -1] .= missing

lat = allowmissing(lat)
lat[lat .== -1] .= missing

lon1 = mapslices(s -> first(skipmissing(s)),lon,dims=2)[:,1]
lat1 = mapslices(s -> first(skipmissing(s)),lat,dims=1)[1,:]


varname = "sbtemper"
sbtemper = nomissing(ds[varname][:,:,:],NaN)

mask = .!(isnan.(sbtemper[:,:,1]) .|| sbtemper[:,:,1] .== 0)

#pcolormesh(mask')

sbtemper[sbtemper .== 0] .= NaN

#clf(); pcolormesh(lon1,lat1,sbtemper[:,:,1]')

mean_sbtemper = mean(sbtemper,dims=3)[:,:,1]
std_sbtemper = std(sbtemper,dims=3)[:,:,1]

ds_Btrc = NCDataset(fnames_Btrc,"r")

@assert ds_Btrc["nav_lon"][:,:] == ds["nav_lon"][:,:]
@assert ds_Btrc["nav_lat"][:,:] == ds["nav_lat"][:,:]



function save(v,varname)
    v = DIVAnd.ufill(v,isfinite.(v))
    interp_fname = joinpath(auxdir,varname * ".nc")
    saveinterp((lon1,lat1),v,
                        (gridlon,gridlat),varname,interp_fname)
end


DOX = nomissing(ds_Btrc["DOX"][:,:,:],NaN)
DOX_min = nomissing(ds_Btrc["DOX_min"][:,:,:],NaN)


low_DOX = mean(DOX_min .< 200,dims=3)[:,:,1]
#pcolor(low_DOX')

mean_DOX = mean(DOX,dims=3)[:,:,1]
std_DOX = std(DOX,dims=3)[:,:,1]


save(mean_sbtemper,"mean_sbtemper")
save(std_sbtemper,"std_sbtemper")

# process DOC

save(mean_DOX,"mean_DOX")
save(std_DOX,"std_DOX")
save(low_DOX,"low_DOX")

# process PAR

PAR = nomissing(ds_Btrc["PAR"][:,:,:],NaN)
mean_PAR = mean(PAR,dims=3)[:,:,1]
std_PAR = std(PAR,dims=3)[:,:,1]

save(mean_PAR,"mean_PAR")
save(std_PAR,"std_PAR")
