from classes import *


def validate_inputs(ice, rt_config, model_config, illumination, impurities):

    """
    Checks all config to make sure all inputs are valid by calling out to:
    validate_snow_algae()
    validate_ilumination()
    validate_ice()


    Args:
    """
    print("\n** Validating model configuration **")
    validate_snow_algae(impurities)
    validate_illumination(illumination)
    validate_ice(ice)

    return 



def validate_snow_algae(impurities):
    """
    Validates snow algae configuration. 
    This includes warning the user that OPs are from literature and not yet field validated.
    Most importantly, checks that the units of the mass absorption coefficient and the cell
    concentration match up (m2/cell -> cells/mL, m2/mg -> ppb)

    Args:
        impurities: a list of Impurity objects
    
    Returns:
        None

    Raises:
        asserts throw exceptions when units mismatch
        otherwise warnings in terminal only
    """
    snw_alg = False
    for i, impurity in enumerate(impurities):
        if ("snw_alg" in impurity.name) and (np.sum(impurity.conc) > 0):
            snw_alg = True
            if impurity.file=="snw_alg_r025um_chla020_chlb025_cara150_carb140.nc":
                if impurity.unit != 0: 
                    raise ValueError("your chosen snow algae file has its MAC in m2/kg, so please express concentration in ppb")
            elif impurity.file=="Data/OP_data/480band/lap/GA_Chevrollier2022_r4.9_L18.8.nc":
                if impurity.unit != 1:
                    raise ValueError("your chosen snow algae file has its MAC in m2/cell, so please express concentration in cells/mL")
            else:
                print("your chosen snow algal OP file is not explicitly validated - please check units carefully")
    
    if snw_alg:
        print("*** WARNING ***")
        print("You have included snow algae as an impurity")
        print("be warned the snow algae optical properties")
        print("are theoretical and yet to be field validated")

        

    print("snow algae OK")

    
    return


def validate_glacier_algae(impurities):
    """
    Validates snow algae configuration. 
    This includes warning the user that OPs are from literature and not yet field validated.
    Most importantly, checks that the units of the mass absorption coefficient and the cell
    concentration match up (m2/cell -> cells/mL, m2/mg -> ppb)

    Args:
        impurities: a list of Impurity objects
    
    Returns:
        None

    Raises:
        asserts throw exceptions when units mismatch
        otherwise warnings in terminal only
    """

    for i, impurity in enumerate(impurities):

        if impurity.file=="GA_Chevrollier2022_r4.9_L18.8.nc":
            if impurity.unit !=0:
                raise ValueError("your chosen glacier algae file has its MAC in m2/cell, so please express concentration in cells/.mL")
        elif impurity.file=="Cook2020_glacier_algae_4_40.nc":
            if impurity.unit !=1:
                raise ValueError("your chosen snow algae file has its MAC in m2/mg, so please express concentration in ppb")
        else:
            print("your chosen glacier algal OP file is not explicitly validated - please check units carefully")

       

    print("glacier algae OK")

    
    return



def validate_illumination(illumination):
    
    if (illumination.solzen > 89) or (illumination.solzen < 1):
        raise ValueError("SZA outside valid range")
    
    if illumination.nbr_wvl != 480:
        raise ValueError("Illumination is has incorrect number of wavelengths")
    
    print("illumination OK")

    return



def validate_ice(ice):
    
    equal_len_fields_ice = [len(ice.dz), len(ice.layer_type), len(ice.cdom), len(ice.rho), len(ice.rds), 
    len(ice.shp), len(ice.water), len(ice.hex_side), len(ice.hex_length), len(ice.ar), ice.nbr_lyr]
    
    if not all(length == equal_len_fields_ice[0] for length in equal_len_fields_ice):
        raise ValueError("variables in ice do not have equal lengths")


    print("ice OK")
    
    return