import numpy as np
import subprocess


def compute_eps2_from_X_and_m( X : float , 
                               m : float , 
                               l_kuhn : float = 1.75e-9 ,
                               N : float = 1100 ,
                               l_omega : float = ( 50e-9 ) * 15 ) -> float:

    return l_kuhn**2 / ( 3/4. * (1-m)*(1+m)*X*l_omega**2 )

def compute_sigma_from_X_and_m( X : float , 
                                m : float , 
                                l_kuhn : float = 1.75e-9 ,
                                N : float = 1100 ,
                                l_omega : float = ( 50e-9 ) * 15 ) -> float:

    return 36 * l_omega**2 / ( ((1-m)/2)**2 * ((1+m)/2)**2 * l_kuhn**2 * X * N**2 )


def main( build_dir : str ):

    m     = 0.5
    X     = 0.3

    eps2  = compute_eps2_from_X_and_m(  X , m )
    sigma = compute_sigma_from_X_and_m( X , m )

    subprocess.run( [ build_dir + "/maDGiCart" , 
                      "--m"     , str(m) ,
                      "--eps2"  , str(eps2) ,
                      "--sigma" , str(sigma) ] )


if __name__ == "__main__":

    build_dir = '../build'

    main( build_dir )