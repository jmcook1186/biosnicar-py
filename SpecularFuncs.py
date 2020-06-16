"""
Author: Joseph Cook, June 2020

contains fresnel equations for calculating specular reflection from upper surface in BioSNICAR_GO

"""

def fresnel(n1, n2, k1, k2, coszen):
        
        """
        function computes fresnel reflection coefficients given the difference in real and imaginary 
        refractive indices between two media (n, k) and the illumination angle, theta
        http://www.oceanopticsbook.info/view/surfaces/the_level_sea_surface
        """
        
        import math
        import numpy as np

        # get zenith angle in degrees from cosine of solar zenith
        theta = np.arccos(coszen) * 180/np.pi

        if theta == 90:
            print("\nIncident angle is 90 degrees")
            Rf = (((n1-1)**2)+k1**2) / (((n2-1)**2)+k2**2)
        
        else:

            theta = math.radians(theta)
            theta_2 = math.asin(1/n2*(math.sin(theta)))
            
            Rf = 0.5*(
                ((math.sin(theta-theta_2)/math.sin(theta+theta_2))**2) +
                ((math.tan(theta-theta_2)/math.tan(theta+theta_2))**2))

        return Rf


