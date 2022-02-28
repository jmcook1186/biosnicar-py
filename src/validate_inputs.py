from classes import *


def validate_inputs(ice, rt_config, model_config, illumination, impurities):

    """
    func checks all config to make sure all inputs are valid
    """
    print("\n** Validating model configuration **")
    check_snow_algae(impurities)
    validate_config_illumination(illumination)
    validate_config_ice(ice)

    return 



def check_snow_algae(impurities):
    """
    If snow algae is present warn user that the OPs are from literature
    """
    snw_alg = False
    for i, impurity in enumerate(impurities):
        if ("snw_alg" in impurity.name) and (np.sum(impurity.conc) > 0):
            snw_alg = True
            if impurity.file=="snw_alg_r025um_chla020_chlb025_cara150_carb140.nc":
                assert(impurity.unit==0, "your chosen snow algae file has its MAC in m2/kg, so please express concentration in ppb")
            elif impurity.file=="Data/OP_data/480band/lap/GA_Chevrollier2022_r4.9_L18.8.nc":
                assert(impurity.unit==1, "your chosen snow algae file has its MAC in m2/cell, so please express concentration in cells/mL")
            else:
                print("your chosen snow algal OP file is not explicitly validated - please check units carefully")
    
    if snw_alg:
        print("*** WARNING ***")
        print("You have included snow algae as an impurity")
        print("be warned the snow algae optical properties")
        print("are theoretical and yet to be field validated")

        

    print("algae OK")

    
    return




def validate_config_illumination(illumination):
    
    assert(illumination.solzen < 89 and illumination.solzen > 1, "SZA outside valid range")
    assert(illumination.nbr_wvl==480, "Illumination is configured for an unusual number of wavelengths")
    
    print("illumination OK")

    return



def validate_config_ice(ice):
    
    equal_len_fields_ice = [len(ice.dz), len(ice.layer_type), len(ice.cdom), len(ice.rho), len(ice.rds), 
    len(ice.shp), len(ice.water), len(ice.hex_side), len(ice.hex_length), len(ice.ar), ice.nbr_lyr]
    
    assert(all(length == equal_len_fields_ice[0] for length in equal_len_fields_ice), "variable length mismatch in ice")    
    assert(ice.ref_idx_re != None, "No ice RI set. Please run ice.calculate_refractive_index()")
    assert(ice.ref_idx_im != None, "No ice RI set. Please run ice.calculate_refractive_index()")
    assert(ice.fl_r_dif_a != None, "No ice RI set. Please run ice.calculate_refractive_index()")
    assert(ice.fl_r_dif_b != None, "No ice RI set. Please run ice.calculate_refractive_index()")
    assert(ice.op_dir !=None, "No ice RI set. Please run ice.calculate_refractive_index()")

    print("ice OK")
    
    return