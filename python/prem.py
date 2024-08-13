def prem(depth):

    """
    isotropic prem model


     given a GLL point, returns super-imposed velocity model values
    """

    HUGEVAL = 1e30
    R_EARTH = 6371000.
    CRUSTAL = False
    SUPPRESS_CRUSTAL_MESH = False
    ONE_CRUST = False

    # GLL point location converted to real
    #xloc = xmesh
    #yloc = ymesh
    #zloc = zmesh

    # get approximate topography elevation at target coordinates
    #distmin = HUGEVAL
    #elevation = 0.0
    #call get_topo_elevation_free_closest(xloc,yloc,elevation,distmin, &
    #                nspec,nglob_dummy,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
    #                num_free_surface_faces,free_surface_ispec,free_surface_ijk)

    # depth in Z-direction
    #if distmin < HUGEVAL:
    #    depth = elevation - zloc
    #else:
    #    depth = - zloc


    # PREM layers (in m)
    R_EARTH_M = 6371000.
    ROCEAN = 6368000.
    RMIDDLE_CRUST = 6356000.
    RMOHO = 6346000.
    R80  = 6291000.
    R220 = 6151000.
    R400 = 5961000. #410
    R600 = 5771000.
    R670 = 5711000. #660
    R771 = 5600000.
    RTOPDDOUBLEPRIME = 3631000.
    RCMB = 3479500.
    RICB = 1217500.

    # compute real physical radius in meters
    r = R_EARTH - depth

    # normalized radius
    x = r / R_EARTH

    # given a normalized radius x, gives the non-dimensionalized density rho,
    # speeds vp and vs, and the quality factors Qkappa and Qmu

    # initializes
    rho = 0.
    vp = 0.
    vs = 0.
    Qmu = 0.
    Qkappa = 0.
    drhodr = 0.

    if (r >= 0. and r <= RICB): # inner core
        drhodr=-2.0*8.8381*x
        rho=13.0885-8.8381*x*x
        vp=11.2622-6.3640*x*x
        vs=3.6678-4.4475*x*x
        Qmu=84.6
        Qkappa=1327.7
    elif (r > RICB and r <= RCMB): # outer core
        drhodr=-1.2638-2.0*3.6426*x-3.0*5.5281*x*x
        rho=12.5815-1.2638*x-3.6426*x*x-5.5281*x*x*x
        vp=11.0487-4.0362*x+4.8023*x*x-13.5732*x*x*x
        vs=0.0
        Qmu=0.0
        Qkappa=57827.0
    elif (r > RCMB and r <= RTOPDDOUBLEPRIME): # D" at the base of the mantle
        drhodr=-6.4761+2.0*5.5283*x-3.0*3.0807*x*x
        rho=7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x
        vp=15.3891-5.3181*x+5.5242*x*x-2.5514*x*x*x
        vs=6.9254+1.4672*x-2.0834*x*x+0.9783*x*x*x
        Qmu=312.0
        Qkappa=57827.0
    elif (r > RTOPDDOUBLEPRIME and r <= R771): # mantle: from top of D" to d670
        drhodr=-6.4761+2.0*5.5283*x-3.0*3.0807*x*x
        rho=7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x
        vp=24.9520-40.4673*x+51.4832*x*x-26.6419*x*x*x
        vs=11.1671-13.7818*x+17.4575*x*x-9.2777*x*x*x
        Qmu=312.0
        Qkappa=57827.0
    elif (r > R771 and r <= R670):
        drhodr=-6.4761+2.0*5.5283*x-3.0*3.0807*x*x
        rho=7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x
        vp=29.2766-23.6027*x+5.5242*x*x-2.5514*x*x*x
        vs=22.3459-17.2473*x-2.0834*x*x+0.9783*x*x*x
        Qmu=312.0
        Qkappa=57827.0
    elif (r > R670 and r <= R600): # mantle: above d670
        drhodr=-1.4836
        rho=5.3197-1.4836*x
        vp=19.0957-9.8672*x
        vs=9.9839-4.9324*x
        Qmu=143.0
        Qkappa=57827.0
    elif (r > R600 and r <= R400):
        drhodr=-8.0298
        rho=11.2494-8.0298*x
        vp=39.7027-32.6166*x
        vs=22.3512-18.5856*x
        Qmu=143.0
        Qkappa=57827.0
    elif (r > R400 and r <= R220):
        drhodr=-3.8045
        rho=7.1089-3.8045*x
        vp=20.3926-12.2569*x
        vs=8.9496-4.4597*x
        Qmu=143.0
        Qkappa=57827.0
    elif (r > R220 and r <= R80):
        drhodr=0.6924
        rho=2.6910+0.6924*x
        vp=4.1875+3.9382*x
        vs=2.1519+2.3481*x
        Qmu=80.0
        Qkappa=57827.0
    else:
        if (CRUSTAL and not SUPPRESS_CRUSTAL_MESH):
            # fill with PREM mantle and later add CRUST2.0
            if (r > R80):
                # density/velocity from mantle just below moho
                drhodr=0.6924
                rho=2.6910+0.6924*x
                vp=4.1875+3.9382*x
                vs=2.1519+2.3481*x
                # shear attenuation for R80 to surface
                Qmu=600.0
                Qkappa=57827.0
        else:
            # use PREM crust
            if (r > R80 and r <= RMOHO):
                drhodr=0.6924
                rho=2.6910+0.6924*x
                vp=4.1875+3.9382*x
                vs=2.1519+2.3481*x
                Qmu=600.0
                Qkappa=57827.0

            elif (SUPPRESS_CRUSTAL_MESH):
                # DK DK extend the Moho up to the surface instead of the crust
                drhodr=0.6924
                rho = 2.6910+0.6924*(RMOHO / R_EARTH)
                vp = 4.1875+3.9382*(RMOHO / R_EARTH)
                vs = 2.1519+2.3481*(RMOHO / R_EARTH)
                Qmu=600.0
                Qkappa=57827.0

            elif (r > RMOHO and r <= RMIDDLE_CRUST):
                drhodr=0.0
                rho=2.9
                vp=6.8
                vs=3.9
                Qmu=600.0
                Qkappa=57827.0

                # same properties everywhere in PREM crust if we decide to define only one layer in the crust
                if (ONE_CRUST):
                    drhodr=0.0
                    rho=2.6
                    vp=5.8
                    vs=3.2
                    Qmu=600.0
                    Qkappa=57827.0


            elif (r > RMIDDLE_CRUST and r <= ROCEAN):
                drhodr=0.0
                rho=2.6
                vp=5.8
                vs=3.2
                Qmu=600.0
                Qkappa=57827.0
            # for density profile for gravity, we do not check that r <= R_EARTH
            elif (r > ROCEAN):
                drhodr=0.0
                rho=2.6
                vp=5.8
                vs=3.2
                Qmu=600.0
                Qkappa=57827.0

    # scales values to SI units ( m/s, kg/m**3)
    #rho_prem = rho*1000.0
    #vp_prem = vp*1000.0
    #vs_prem = vs*1000.0

    qmu_atten = Qmu
    qkappa_atten = Qkappa

    return rho, vp, vs, qmu_atten, qkappa_atten
