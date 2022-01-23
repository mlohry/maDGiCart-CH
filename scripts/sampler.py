import numpy as np
import subprocess
import csv
import sys , os


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


def compute_eps2_and_sigma_from_X_and_m_3D( X : float , 
                                            m : float , 
                                            l_kuhn : float = 1.75e-9 ,
                                            N : float = 1100 ,
                                            l_omega : float = ( 50e-9 ) * 15 ) -> float:

    l_omega_2D = l_omega * 1e9    # domain side length of 2D square (in nm)

    l_omega_3D = ( l_omega_2D * l_omega_2D * 1 )**(1/3) # 3D cube length that gives same volume as 2D square area
    
    l_omega_3D *= 1e-9 # scaling to get nm
    
    return compute_eps2_from_X_and_m(  X , m , l_kuhn , N , l_omega_3D ) , compute_sigma_from_X_and_m( X , m , l_kuhn , N , l_omega_3D )


def main( outdir : str , builddir : str , nsamples : int ):

    m     = np.random.uniform( -0.7  , 0.7 , nsamples )
    X     = np.random.uniform( 0.055 , 0.5 , nsamples )

    for i,(mi,Xi) in enumerate( zip(m,X) ):

        workdir = outdir + '/run_' + str(i+1)
        subprocess.run( [ 'mkdir' , workdir ] )
        os.chdir( workdir )

        eps2  = compute_eps2_from_X_and_m(  Xi , mi )
        sigma = compute_sigma_from_X_and_m( Xi , mi )

        subprocess.run( [ builddir + "/maDGiCart" ,
                        "--m"     , str(mi) ,
                        "--eps2"  , str(eps2) ,
                        "--sigma" , str(sigma) ] )

        # Save mc parameters
        w = csv.writer( open( workdir + "/params_" + str(i+1) + ".csv" , "w" ) )
        w.writerow( [ 'm' , mi ] )
        w.writerow( [ 'X' , Xi ] )



if __name__ == "__main__":

    outdir    = sys.argv[1]
    builddir  = sys.argv[2]
    nsamples  = int( sys.argv[3] )

    main( outdir , builddir , nsamples )
