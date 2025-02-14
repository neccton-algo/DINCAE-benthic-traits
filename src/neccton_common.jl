# Common functions to all steps
# in the preparing, training and validation of the reconstruction

using NCDatasets
using Interpolations

# Resolution
dlon = dlat = 1/16

# lon/lat grid
gridlon = 27.4:dlon:32
gridlat = 41.8:dlat:46.8

# Directories
basedir = expanduser("~/Reconstruct_Points/Datasets")
moddir = expanduser("~/Reconstruct_Points/ModelOut")
moddirNew = expanduser("~/Reconstruct_Points/ModelNew")
figdir = expanduser("~/Reconstruct_Points/Datasets/Figures")

datadir = joinpath(basedir,"BenthicData")
auxdir = joinpath(basedir,"Auxdata_$(1/dlon)")

# Files
split_fname = joinpath(basedir,"split.nc")

station_fname = joinpath(datadir,"Stations","stations.csv")
CWM_response_fname = joinpath(datadir,"CWM_SxT","cwm.csv")
env_matrix_fname = joinpath(datadir,"Environnement","matrice environnement.txt")


function neccton_load_mask()

    #fnames0 = sort(glob("Cl*_1d_*_*_btrc_T_*-*.nc",moddir))
    fnames = sort(glob("Cl*_1d_*_*_grid_T_*-*.nc",moddir))
    NCDataset(fnames[1]) do ds_T
        sbtemper = nomissing(ds_T["sbtemper"][:],NaN)
        mask = .!(isnan.(sbtemper[:,:,1]) .|| sbtemper[:,:,1] .== 0)
        return mask
    end
end


"""
    lon,lat = neccton_load_coord(fname)

Load longitude and latitude from NEMO file `fname`.
"""
function neccton_load_coord(fname)
        ds = NCDataset(fname,"r")
        lon = ds["nav_lon"][:,:]
        lat = ds["nav_lat"][:,:]

        lon = allowmissing(lon)
        lon[lon .== -1] .= missing

        lat = allowmissing(lat)
        lat[lat .== -1] .= missing

        lon1 = mapslices(s -> first(skipmissing(s)),lon,dims=2)[:,1]
        lat1 = mapslices(s -> first(skipmissing(s)),lat,dims=1)[1,:]

        return lon1,lat1
    end


"""
    df = df_load(station_fname,CWM_response_fname)

Load dataframe from file `station_fname` and `CWM_response_fname` and join
on station name.
"""
function df_load(station_fname,CWM_response_fname)
    station = CSV.read(station_fname,DataFrame);
    resp = CSV.read(CWM_response_fname,DataFrame,decimal=',');
    resp = rename(resp, :Column1 => :sta)
    df = innerjoin(resp, station, on = :sta)
    return df, resp, station
end


"""
    saveinterp((lon,lat),field,(gridlon,gridlat),varname,interp_fname)

Interpolate the field `field` defined over the grid `(lon,lat)`
on the grid `(gridlon,gridlat)`. The result is saved in the NetCDF files
`interp_fname` under the same `varname`.
`lon`, `lat`, `gridlon`, `gridlat` are all vectors.
"""
function saveinterp((lon,lat),SS2,(gridlon,gridlat),varname,interp_fname)
    @info "interpolate"
    itp = interpolate((lon,lat), SS2, Gridded(Linear()));
    SSi = itp(gridlon,gridlat);

    @show extrema(SSi)

    ds = Dataset(interp_fname,"c")
    # Dimensions

    ds.dim["lon"] = length(gridlon)
    ds.dim["lat"] = length(gridlat)

    # Declare variables

    nclon = defVar(ds,"lon", Float64, ("lon",))
    nclon.attrib["units"] = "degrees_east"
    nclon.attrib["standard_name"] = "longitude"
    nclon.attrib["long_name"] = "longitude"

    nclat = defVar(ds,"lat", Float64, ("lat",))
    nclat.attrib["units"] = "degrees_north"
    nclat.attrib["standard_name"] = "latitude"
    nclat.attrib["long_name"] = "latitude"

    ncvar = defVar(ds,lowercase(varname), Float32, ("lon", "lat"))
    ncvar.attrib["_FillValue"] = Float32(9.96921e36)
    ncvar.attrib["missing_value"] = Float32(9.96921e36)
    ncvar.attrib["long_name"] = varname


    # Define variables

    nclon[:] = gridlon
    nclat[:] = gridlat
    ncvar[:,:] = SSi

    close(ds)
end
