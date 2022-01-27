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
    
    l_omega_3D *= 10 # annoying empirical fudge factor, need to think about this more... 

    l_omega_3D *= 1e-9 # scaling to get nm
    
    return compute_eps2_from_X_and_m(  X , m , l_kuhn , N , l_omega_3D ) , compute_sigma_from_X_and_m( X , m , l_kuhn , N , l_omega_3D )


def main_2d_single_run( outdir : str , builddir : str ):

    mi     = np.random.uniform( -0.7  , 0.7 )
    Xi     = np.random.uniform( 0.055 , 0.5 )

    os.chdir( outdir )

    with open( outdir + "/params.csv" , 'w' ) as myfile:
        w = csv.writer( myfile )
        w.writerow( [ 'm' , mi ] )
        w.writerow( [ 'X' , Xi ] )
        myfile.flush()

    eps2  = compute_eps2_from_X_and_m(  Xi , mi )
    sigma = compute_sigma_from_X_and_m( Xi , mi )

    subprocess.run( [ builddir + "/maDGiCart" ,
                    "--m"     , str(mi) ,
                    "--eps2"  , str(eps2) ,
                    "--sigma" , str(sigma) ] )


def main_3d_single_run( outdir : str , builddir : str ):

    mi     = np.random.uniform( -0.7  , 0.7 )
    Xi     = np.random.uniform( 0.055 , 0.5 )
    
    os.chdir( outdir )
    
    with open( outdir + "/params.csv" , 'w' ) as myfile:
        w = csv.writer( myfile )
        w.writerow( [ 'm' , mi ] )
        w.writerow( [ 'X' , Xi ] )
        myfile.flush()

    eps2 , sigma = compute_eps2_and_sigma_from_X_and_m_3D( Xi , mi )

    subprocess.run( [ builddir + "/maDGiCart" ,
                        "--dimension" , "3" ,
                        "--max_time_steps"  , "100000" , 
                        "--time_integrator" , "ode23" ,
                        "--converged_rel_tol" , "1e-6" ,
                        "--m"     , str(mi) ,
                        "--eps2"  , str(eps2) ,
                        "--sigma" , str(sigma) ] )


def main_3d_single_run_from_file( outdir : str , builddir : str , m : float , X : float , initial_condition_file : str ):

    os.chdir( outdir )
    
    with open( outdir + "/params.csv" , 'w' ) as myfile:
        w = csv.writer( myfile )
        w.writerow( [ 'm' , m ] )
        w.writerow( [ 'X' , X ] )
        myfile.flush()

    eps2 , sigma = compute_eps2_and_sigma_from_X_and_m_3D( X , m )

    subprocess.run( [ builddir + "/maDGiCart" ,
                        "--dimension" , "3" ,
                        "--max_time_steps"  , "10000" , 
                        "--time_integrator" , "ode23" ,
                        "--converged_rel_tol" , "1e-4" ,
                        "--initial_condition_file" , initial_condition_file ,
                        "--m"     , str(m) ,
                        "--eps2"  , str(eps2) ,
                        "--sigma" , str(sigma) ] )



if __name__ == "__main__":

    #outdir    = sys.argv[1]
    #builddir  = sys.argv[2]

    #main_3d_single_run( outdir , builddir )
    
    outdir   = '/home/adegennaro/Projects/appmath/madgicartCH-postproc/data/stacked_m0p4_X0p15'
    builddir = '/home/adegennaro/Projects/appmath/maDGiCart-CH-build'

    main_3d_single_run_from_file( outdir , builddir , 0.4 , 0.15305960080040307 , 
        '/home/adegennaro/Projects/appmath/madgicartCH-postproc/data/stacked_m0p4_X0p15/c_stack_cube_m0p4_ascii.vts' )

