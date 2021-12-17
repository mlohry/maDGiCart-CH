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


def main( topdir : str , builddir : str , outdir : str , nsamples : int ):

    m     = np.random.uniform( -0.7  , 0.7 , nsamples )
    X     = np.random.uniform( 0.055 , 0.5 , nsamples )

    for i,(mi,Xi) in enumerate( zip(m,X) ):

        workdir = topdir + '/' + outdir + '/run_' + str(i+1)
        subprocess.run( [ 'mkdir' , workdir ] )
        os.chdir( workdir )

        eps2  = compute_eps2_from_X_and_m(  Xi , mi )
        sigma = compute_sigma_from_X_and_m( Xi , mi )

        subprocess.run( [ topdir + '/' + builddir + "/maDGiCart" ,
                        "--interactive_plots" , 'false' ,  
                        "--m"     , str(mi) ,
                        "--eps2"  , str(eps2) ,
                        "--sigma" , str(sigma) ] )

        # Save mc parameters
        w = csv.writer( open( workdir + "/params_" + str(i+1) + ".csv" , "w" ) )
        w.writerow( [ 'm' , mi ] )
        w.writerow( [ 'X' , Xi ] )



if __name__ == "__main__":

    topdir    = sys.argv[1]
    builddir = sys.argv[2]
    outdir    = sys.argv[3]
    nsamples  = int( sys.argv[4] )

    main( topdir , builddir , outdir , nsamples )
