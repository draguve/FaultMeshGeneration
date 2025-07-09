import math

import numpy as np
import typer
from netCDF4 import Dataset


def writeNetcdf4Paraview(sname, x, y, z, aName, aData):
    "create a netcdf file readable by paraview (but not by ASAGI)"
    fname = sname + "_paraview.nc"
    print("writing " + fname)
    # Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]
    rootgrp = Dataset(fname, "w", format="NETCDF4")
    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)
    rootgrp.createDimension("w", nz)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    vz = rootgrp.createVariable("w", "f4", ("w",))
    vz[:] = z
    for i in range(len(aName)):
        vTd = rootgrp.createVariable(aName[i], "f4", ("u", "v", "w"))
        vTd[:, :, :] = aData[i][:, :, :]
    rootgrp.close()


def writeNetcdf4SeisSol(sname, x, y, z, aName, aData):
    # "create a netcdf file readable by ASAGI (but not by paraview)"
    # creating the file for SeisSol
    fname = sname + "_ASAGI.nc"
    print("writing " + fname)
    # Creating the netcdf file
    nx = x.shape[0]
    ny = y.shape[0]
    nz = z.shape[0]
    rootgrp = Dataset(fname, "w", format="NETCDF4")

    rootgrp.createDimension("u", nx)
    rootgrp.createDimension("v", ny)
    rootgrp.createDimension("w", nz)

    vx = rootgrp.createVariable("u", "f4", ("u",))
    vx[:] = x
    vy = rootgrp.createVariable("v", "f4", ("v",))
    vy[:] = y
    vz = rootgrp.createVariable("w", "f4", ("w",))
    vz[:] = z
    ldata4 = [(name, "f4") for name in aName]
    ldata8 = [(name, "f8") for name in aName]
    mattype4 = np.dtype(ldata4)
    mattype8 = np.dtype(ldata8)
    mat_t = rootgrp.createCompoundType(mattype4, "material")

    # this transform the 4 D array into an array of tuples
    arr = np.stack([aData[i] for i in range(len(aName))], axis=3)
    newarr = arr.view(dtype=mattype8)
    newarr = newarr.reshape(newarr.shape[:-1])
    mat = rootgrp.createVariable("data", mat_t, ("u", "v", "w"))
    mat[:] = newarr
    rootgrp.close()

def main():
    x = np.linspace(-50000, 50000, int((50000 + 50000) / 100) + 1)  # for testing if correct
    y = np.linspace(0, 10000, int(10000 / 100))
    z = np.linspace(-50000, 50000, int((50000 + 50000) / 100))

    # Generate a Gaussian shaped slip distribution for illustration purposes
    xg, yg, zg = np.meshgrid(x, y, z, indexing='ij')
    Rx = np.zeros(xg.shape)
    Rz = np.zeros(xg.shape)

    loc_temp = xg < -9800
    Rx[loc_temp] = (-xg[loc_temp] - 9800) / 10000
    loc_temp = xg > 1100
    Rx[loc_temp] = (xg[loc_temp] - 1100) / 10000

    loc_temp = zg < -8000
    Rz[loc_temp] = (-zg[loc_temp] - 8000) / 1000
    loc_temp = zg > -2300
    Rz[loc_temp] = (zg[loc_temp] + 2300) / 10000

    Rt = np.minimum(1, np.sqrt(np.power(Rx, 2) + np.power(Rz, 2)))
    s_xy0 = 30000000.0 * (1.0 - Rt)

    radius = np.sqrt(np.power(xg + 6000, 2) + np.power(zg + 6000, 2))
    pi = 4 * math.atan(1.0)

    s_xy1 = np.zeros(xg.shape)
    loc_temp = radius <= 550
    s_xy1[loc_temp] = 3150000
    loc_temp = 1575000. * (1.0 + np.cos(pi * (radius - 550.0) / 250.0))

    sheer_stress = s_xy0 + s_xy1

    prefix = "test"
    ldataName = ["sheer_stress"]
    lgridded_myData = [sheer_stress]

    writeNetcdf4Paraview(prefix, x, y, z, ldataName, lgridded_myData)
    writeNetcdf4SeisSol(prefix, x, y, z, ldataName, lgridded_myData)


if __name__ == '__main__':
    typer.run(main)