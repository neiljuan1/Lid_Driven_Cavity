import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import sparse
from scipy.sparse.linalg import gmres
from scipy.sparse.linalg import lgmres
from scipy.sparse.linalg import cgs


def main():
    import time 
    import post_processor
    import numpy as np
    import math

    start_time = time.time()

    # Ghia et al. Conditions
    U = 1 # lid velocity
    L = 1 # domain length
    Re = 100
    rho = 1
    mu = rho*L/Re
    t = 1
    #npts = 5
    nx = 257
    ny = 257

    dy = dx = L/(nx-1)

    # domain for postprocessing
    x = np.zeros(nx)
    y = np.zeros(ny)

    for i in range(1,len(x)):
        x[i]=x[i-1] + dx
        y[i]=y[i-1] + dx


        # Error is the momentum residual at all velocity nodes 
    # and monitor the sum of absolute values of the residuals.
    error = 1
    errcon = 1e-6
    u_res_vec, v_res_vec, p_res_vec = [],[],[]
    u_res = np.zeros(shape = (nx,ny+1))
    v_res = np.zeros(shape = (nx+1,ny))
    p_res = np.zeros(shape = (nx+1,ny+1))


    alpha_u = 0.8 # URF
    alpha_p = 0.8
    alpha_v = alpha_u
    maxiter = 50000
    niter = 1

    # Collocated Solution Fields
    u_sol = np.zeros(shape = (nx,ny))
    v_sol = np.zeros(shape = (nx,ny)) # final u nodes i = 1,2,3,4
    p_sol = np.zeros(shape = (nx,ny)) # final p nodes A,B,C,D,E


    # Final
    u = np.zeros(shape = (nx,ny+1))
    v = np.zeros(shape = (nx+1,ny)) # final u nodes i = 1,2,3,4
    p = np.zeros(shape = (nx+1,ny+1)) # final p nodes A,B,C,D,E


    # Guessed Fields East is [0,:], West is [-1,:], South is [:,-1], North is [:,0]
    u_star = np.zeros(shape = (nx,ny+1))
    v_star = np.zeros(shape = (nx+1,ny)) # final u nodes i = 1,2,3,4
    p_star = np.zeros(shape = (nx+1,ny+1)) # final p nodes A,B,C,D,E

    u_old = np.zeros(shape = (nx,ny+1))
    v_old = np.zeros(shape = (nx+1,ny)) # final u nodes i = 1,2,3,4
    p_old = np.zeros(shape = (nx+1,ny+1)) # final p nodes A,B,C,D,E


    # Discreiztization Variables 
    # U velocity Field Coefficients
    ua_E = np.zeros(shape = (nx,ny+1))
    ua_W = np.zeros(shape = (nx,ny+1))
    ua_N = np.zeros(shape = (nx,ny+1))
    ua_S = np.zeros(shape = (nx,ny+1))
    ua_P = np.zeros(shape = (nx,ny+1))

    uS_u = np.zeros(shape = (nx,ny+1))
    ud = np.zeros(shape = (nx,ny+1))

    # V velocity Field Coefficients
    va_E = np.zeros(shape = (nx+1,ny))
    va_W = np.zeros(shape = (nx+1,ny))
    va_N = np.zeros(shape = (nx+1,ny))
    va_S = np.zeros(shape = (nx+1,ny))
    va_P = np.zeros(shape = (nx+1,ny))

    vS_u = np.zeros(shape = (nx+1,ny))
    vd = np.zeros(shape = (nx+1,ny))

    # Pressure Field Coefficients
    pa_E = np.zeros(shape = (nx+1,ny+1))
    pa_W = np.zeros(shape = (nx+1,ny+1))
    pa_N = np.zeros(shape = (nx+1,ny+1))
    pa_S = np.zeros(shape = (nx+1,ny+1))
    pa_P = np.zeros(shape = (nx+1,ny+1)) 
    pbp = np.zeros(shape = (nx+1,ny+1)) # b prime
    pp = np.zeros(shape = (nx+1,ny+1)) # p prime field


    # Boundary Conditions
    # U
    # East
    u[0,:] = 0
    v[0,:] = -v[1,:] #0
    p[0,:] = p[1,:]
    # West
    u[-1,:] = 0
    v[-1,:] = -v[-2,:] #0
    p[-1,:] = p[-2,:]
    # North 
    u[:,0] = U # lid
    v[:,0] = 0 # lid
    p[:,0] = p[:,1] # lid
    # South
    u[:,-1] = -u[:,-2]
    v[:,-1] = 0
    p[:,-1] = p[:,-2]


    # Initial Guesses
    p[:,:] = 0#1
    u_star = u.copy()
    v_star = v.copy()

    # Save old solution 
    u_old = u.copy()
    v_old = v.copy()
    p_old = p.copy()


    # original
    while (error>errcon and niter<maxiter):
        # Step 1 Solve the discretized momentum equations. [row, col] = [E -> W, N -> S] OBTAIN u* v* 
        # (a_P*u_star_P = a_W*u_star_W*+a_E*u_star_E + Su)
        # Su = (P_W - P_E)*A_w

        # Use Central Differencing Scheme. Need a coefficients

        # Interior u nodes
        # Convective F's F = rho*u_face*A_face
        # velocity at faces:
        for i in range(1,nx-1):
            for j in range(1,ny): 
                u_star_e = (u_old[i-1,j] + u_old[i,j])/2
                u_star_w = (u_old[i+1,j] + u_old[i,j])/2
                v_star_n = (v_old[i,j-1] + v_old[i+1,j-1])/2
                v_star_s = (v_old[i,j] + v_old[i+1,j])/2

                F_w = rho * u_star_w * dy*t
                F_e = rho * u_star_e * dy*t
                F_n = rho * v_star_n * dx*t
                F_s = rho * v_star_s * dx*t

                # Diffusive conductance
                D_w = mu / dx * dy*t
                D_e = mu / dx * dy*t
                D_n = mu / dy * dy*t
                D_s = mu / dy * dx*t

                # Hybrid
                ua_W[i,j] = max(F_w, D_w + F_w/2, 0)
                ua_E[i,j] = max(-F_e, D_e - F_e/2, 0)
                ua_N[i,j] = max(-F_n, D_n - F_n/2, 0)
                ua_S[i,j] = max(F_s, D_s + F_s/2, 0)

                ua_P[i,j] = (ua_E[i,j] + ua_W[i,j] + ua_N[i,j] + ua_S[i,j])

                uS_u[i,j] = (p_old[i+1,j]-p_old[i,j])*dy*t

                ud[i,j] = dy*t/ua_P[i,j]


        # Solve velocity field from momentum equationfor i in range(1,nx-1):
        for i in range(1,nx-1):
            for j in range(1,ny): 
                u_star[i,j] = (ua_W[i,j]*u_old[i+1,j] + ua_E[i,j]*u_old[i-1,j] + ua_N[i,j]*u_old[i,j-1] + ua_S[i,j]*u_old[i,j+1] + \
                            uS_u[i,j])/ua_P[i,j]
        
        
        # Interior v nodes
        # Convective F's F = rho*u_face*A_face
        # velocity at faces:
        for i in range(1,nx):
            for j in range(1,ny-1):
                u_star_e = (u_old[i-1,j] + u_old[i-1,j+1])/2
                u_star_w = (u_old[i,j] + u_old[i,j+1])/2
                v_star_n = (v_old[i,j-1] + v_old[i,j])/2
                v_star_s = (v_old[i,j+1] + v_old[i,j])/2

                F_w = rho * u_star_w * dy*t
                F_e = rho * u_star_e * dy*t
                F_n = rho * v_star_n * dx*t
                F_s = rho * v_star_s * dx*t

                # Diffusive conductance
                D_w = mu / dx * dy*t
                D_e = mu / dx * dy*t
                D_n = mu / dy * dy*t
                D_s = mu / dy * dx*t

                # Central Difference
                #va_W[i,j] = D_w + F_w/2
                #va_E[i,j] = D_e - F_e/2
                #va_N[i,j] = D_n - F_n/2
                #va_S[i,j] = D_s + F_s/2

                # Hybrid
                va_W[i,j] = max(F_w, D_w + F_w/2, 0)
                va_E[i,j] = max(-F_e, D_e - F_e/2, 0)
                va_N[i,j] = max(-F_n, D_n - F_n/2, 0)
                va_S[i,j] = max(F_s, D_s + F_s/2, 0)

                va_P[i,j] = (va_E[i,j] + va_W[i,j] + va_N[i,j] + va_S[i,j])

                vS_u[i,j] = (p_old[i,j+1]-p_old[i,j])*dy*t 

                vd[i,j] = dy*t/va_P[i,j]

        # Solve velocity field from momentum equationfor i in range(1,nx-1):
        for i in range(1,nx):
            for j in range(1,ny-1): 
                v_star[i,j] = (va_W[i,j]*v_old[i+1,j] + va_E[i,j]*v_old[i-1,j] + va_N[i,j]*v_old[i,j-1] + va_S[i,j]*v_old[i,j+1] + \
                            vS_u[i,j])/va_P[i,j]      
        
        # Step 2: Solve the Pressure Correction Equation OBTAIN P'
        # Find Coefficients using the continuity equation

        for i in range(1,nx):
            for j in range(1,ny):
                pa_W[i,j] = rho * ud[i,j] * dy*t
                pa_E[i,j] = rho * ud[i-1,j] * dy*t
                pa_N[i,j] = rho * vd[i,j-1] * dx*t
                pa_S[i,j] = rho * vd[i,j] * dx*t

                pa_P[i,j] = pa_W[i,j] + pa_E[i,j] + pa_N[i,j] + pa_S[i,j]


                pF_w_star = rho * u_star[i,j] * dy*t
                pF_e_star = rho * u_star[i-1,j] * dy*t

                pF_n_star = rho * v_star[i,j-1] * dx*t
                pF_s_star = rho * v_star[i,j] * dx*t

                # Flux in - flux out. Fw - Fe + Fs- Fn
                pbp[i,j] = (pF_w_star - pF_e_star) + (pF_s_star - pF_n_star) # b prime

        pp[:,:] = 0
        for i in range(1,nx):
            for j in range(1,ny):
                pp[i,j] = (pa_W[i,j]*pp[i+1,j] + pa_E[i,j]*pp[i-1,j] + pa_S[i,j]*pp[i,j+1] + pa_N[i,j]*pp[i,j-1] \
                        + pbp[i,j])/pa_P[i,j]
                
                
                
        # Step 3: Correct pressure and velocities ! P_S - P_P and P_W - P_P

        p = p_old + pp * alpha_p
        for i in range(1,nx-1):
            for j in range(1,ny): 
                u[i,j] = u_star[i,j] + ud[i,j]*(pp[i+1,j] - pp[i,j])

        for i in range(1,nx):
            for j in range(1,ny-1): 
                v[i,j] = v_star[i,j] + vd[i,j]*(pp[i,j+1] - pp[i,j]) 

        # Under-relax for the next solution field
        
        u = (1-alpha_u)*u_old + alpha_u*u
        v = (1-alpha_v)*v_old + alpha_v*v

        # Boundary Conditions
        # U
        # East
        u[0,:] = 0
        v[0,:] = -v[1,:]
        p[0,:] = p[1,:]
        # West
        u[-1,:] = 0
        v[-1,:] = -v[-2,:]
        p[-1,:] = p[-2,:]
        # North 
        u[:,0] =  U#2 - u[:,1]# lid
        v[:,0] = 0 # lid
        p[:,0] = p[:,1] # lid
        # South
        u[:,-1] = -u[:,-2]#0
        v[:,-1] = 0
        p[:,-1] = p[:,-2]

        # Error Analysis
        u_res = abs(u - u_old)
        v_res = abs(v - v_old)
        p_res = abs(p - p_old)

        max_u_res = max(map(max,u_res))
        max_v_res = max(map(max,v_res))
        max_p_res = max(map(max,p_res))
        u_res_vec.append(math.log(max_u_res))
        v_res_vec.append(math.log(max_v_res))
        p_res_vec.append(math.log(max_p_res))
        error = max(max_u_res,max_v_res,max_p_res)

        # Next iteration - Set guessed values as the past iteration solution
        u_old = u.copy()
        v_old = v.copy()
        p_old = p.copy()

        niter+=1
        print("_______iter______{}_____error__________________{}".format(niter,error))

    for i in range(0,nx):
        for j in range(0,ny):
            u_sol[i,j] = (u[i,j] + u[i,j+1])/2
            v_sol[i,j] = (v[i,j] + v[i+1,j])/2
            p_sol[i,j] = (p[i,j] + p[i,j+1] + p[i+1,j] + p[i+1,j+1])/4

    u_sol_new = u_sol.copy()
    v_sol_new = v_sol.copy()
    p_sol_new = p_sol.copy()

    for i in range(nx):
        for j in range(ny):
            u_sol_new[i,j] = u_sol[ny-j-1,i]
            v_sol_new[i,j] = v_sol[ny-j-1,i]
            p_sol_new[i,j] = p_sol[ny-j-1,i]
        

    post_processor.contour_plot(x,y,nx,ny,u_sol_new,v_sol_new,p_sol_new)
    post_processor.horizontal(x,y,nx,ny,u_sol_new,v_sol_new,p_sol_new)
    post_processor.vertical(x,y,nx,ny,u_sol_new,v_sol_new,p_sol_new)

    print('time elapsed: ',time.time()-start_time)
    plt.plot(u_res_vec,'r',label = "U Residual")
    plt.plot(v_res_vec,'b',label = "V Residual")
    plt.plot(p_res_vec,'g',label = "P Residual")
    plt.xlabel('Iterations')
    plt.ylabel('Residual Error')
    plt.legend()
    plt.show()
    
    

if __name__=='__main__' :
    main()