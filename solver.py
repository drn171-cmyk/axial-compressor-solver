# AE3007 Compressor Aerothermal Analysis
# LPC (Low-Pressure Compressor) + HPC (High-Pressure Compressor)
# 5 radial positions per stage: 0=hub, 1=lower-mid, 2=mean, 3=upper-mid, 4=tip

import numpy as np

# CONSTANTS & INPUTS
massFlowRate=110 # kg/s
Cp=1005 # J/kgK
gamma=1.4 
R=287 # J/kgK
blockage=0.7 # annulus blockage factor

Vx=277.34 # m/s
pInTot=50535.6 # Pa
TempInTot=624.45 # K
RpmLow=8700 # RPM
RpmHigh=16123 # RPM

stageNumLow=6
stageNumHigh=8

rHub=0.2 #m

# MANUAL BETA ANGLES (Degrees)
beta2_mean_deg_LPC=48
beta2_mean_LPC=np.deg2rad(beta2_mean_deg_LPC)

beta2_mean_deg_HPC=58
beta2_mean_HPC=np.deg2rad(beta2_mean_deg_HPC)

# INITIAL CALCULATIONS
omegaLow=2*np.pi*RpmLow/60
omegaHigh=2*np.pi*RpmHigh/60

tol=1e-6 # tolerance
err=1 # Error

TempIn=TempInTot-Vx**2/(2*Cp)
pIn=pInTot*(TempIn/TempInTot)**(gamma/(gamma-1))
rhoIn=pIn/(R*TempIn)


# CORE ANALYSIS FUNCTION
def analyse_spool(stageNum, omega, pTot0, TTot0, beta2_mean_spool, label):
    nStat=stageNum*4
    nPlane=stageNum*2+1

    Vel=np.zeros((5,nStat))
    U=np.zeros((5,nStat))
    W=np.zeros((5,nStat))
    VelTheta=np.zeros((5,nStat))
    alpha=np.zeros((5,nStat))
    beta=np.zeros((5,nStat))

    pTot=np.zeros((5,nPlane))
    TTot=np.zeros((5,nPlane))
    pStat=np.zeros((5,nPlane))
    TStat=np.zeros((5,nPlane))
    rho=np.zeros((5,nPlane))

    rTip=np.zeros(stageNum)
    rMean=np.zeros(stageNum)

    TStat0=TTot0-Vx**2/(2*Cp)
    pStat0=pTot0*(TStat0/TTot0)**(gamma/(gamma-1))
    rho0=pStat0/(R*TStat0)

    pTot[:,0]=pTot0
    TTot[:,0]=TTot0
    pStat[:,0]=pStat0
    TStat[:,0]=TStat0
    rho[:,0]=rho0

    for s in range(stageNum):
        plane_in=2*s
        plane_rot=2*s+1
        plane_out=2*s+2

        stat_ri=s*4
        stat_ro=s*4+1
        stat_si=s*4+2
        stat_so=s*4+3

        rho_in=rho[2,plane_in]
        rTip[s]=np.sqrt(massFlowRate/(blockage*np.pi*rho_in*Vx)+rHub**2)
        rMean[s]=rHub+(rTip[s]-rHub)/2.0

        r_m=rMean[s]
        U_m=omega*r_m
        W_m=Vx/np.cos(beta2_mean_spool)
        Vtheta_m=U_m-W_m*np.sin(beta2_mean_spool)

        for i in range(5):
            r_i=rHub+(rTip[s]-rHub)*i/4.0
            U_i=omega*r_i

            # Rotor inlet
            Vel[i,stat_ri]=Vx
            U[i,stat_ri]=U_i
            W[i,stat_ri]=np.sqrt(Vx**2+U_i**2)
            VelTheta[i,stat_ri]=0.0
            beta[i,stat_ri]=np.arctan(U_i/Vx)
            alpha[i,stat_ri]=0.0

            # Rotor exit
            Vtheta_i=Vtheta_m*r_m/r_i
            beta2_i=np.arctan((U_i-Vtheta_i)/Vx)

            U[i,stat_ro]=U_i
            beta[i,stat_ro]=beta2_i
            W[i,stat_ro]=Vx/np.cos(beta2_i)
            VelTheta[i,stat_ro]=Vtheta_i
            Vel[i,stat_ro]=np.sqrt(Vx**2+Vtheta_i**2)
            alpha[i,stat_ro]=np.arctan(Vtheta_i/Vx)

            # Rotor thermodynamics
            dH=U_i*Vtheta_i
            T_tot_in=TTot[i,plane_in]
            p_tot_in=pTot[i,plane_in]

            TTot[i,plane_rot]=T_tot_in+dH/Cp
            pTot[i,plane_rot]=p_tot_in*(TTot[i,plane_rot]/T_tot_in)**(gamma/(gamma-1))
            TStat[i,plane_rot]=TTot[i,plane_rot]-Vel[i,stat_ro]**2/(2*Cp)
            pStat[i,plane_rot]=pTot[i,plane_rot]*(TStat[i,plane_rot]/TTot[i,plane_rot])**(gamma/(gamma-1))
            rho[i,plane_rot]=pStat[i,plane_rot]/(R*TStat[i,plane_rot])

            # Stator inlet
            Vel[i,stat_si]=Vel[i,stat_ro]
            alpha[i,stat_si]=alpha[i,stat_ro]
            VelTheta[i,stat_si]=VelTheta[i,stat_ro]

            # Stator exit
            alpha[i,stat_so]=0.0
            Vel[i,stat_so]=Vx
            VelTheta[i,stat_so]=0.0

            # Stator thermodynamics
            TTot[i,plane_out]=TTot[i,plane_rot]
            pTot[i,plane_out]=pTot[i,plane_rot]
            TStat[i,plane_out]=TTot[i,plane_out]-Vx**2/(2*Cp)
            pStat[i,plane_out]=pTot[i,plane_out]*(TStat[i,plane_out]/TTot[i,plane_out])**(gamma/(gamma-1))
            rho[i,plane_out]=pStat[i,plane_out]/(R*TStat[i,plane_out])

    return dict(Vel=Vel, U=U, W=W, VelTheta=VelTheta, alpha=alpha, beta=beta, pTot=pTot, TTot=TTot, pStat=pStat, TStat=TStat, rho=rho, rTip=rTip, rMean=rMean, stageNum=stageNum, omega=omega, label=label)


# OUTPUT & CAD EXPORT FUNCTIONS
def print_geometry_for_cad(d, stage_index=0):
    label=d['label']
    print(f"\n=== {label} STAGE {stage_index+1} BLADE GEOMETRY FOR CAD ===")
    print(f"{'Section':<10} | {'r (m)':<8} | {'Beta1 (deg)':<12} | {'Beta2 (deg)':<12} | {'Alpha2 (deg)':<12} | {'Camber (deg)':<12}")
    print("-" * 78)

    stat_ri=stage_index*4
    stat_ro=stage_index*4+1
    sections=["Hub", "25%", "Mean", "75%", "Tip"]
    
    for i in range(5):
        r_val=rHub+(d['rTip'][stage_index]-rHub)*i/4.0
        b1=np.degrees(d['beta'][i,stat_ri])
        b2=np.degrees(d['beta'][i,stat_ro])
        a2=np.degrees(d['alpha'][i,stat_ro])
        camber=b1-b2
        print(f"{sections[i]:<10} | {r_val:<8.4f} | {b1:<12.2f} | {b2:<12.2f} | {a2:<12.2f} | {camber:<12.2f}")


def print_stator_geometry_for_cad(d, stage_index=0):
    label=d['label']
    print(f"\n=== {label} STAGE {stage_index+1} STATOR GEOMETRY FOR CAD ===")
    print(f"{'Section':<10} | {'r (m)':<8} | {'Alpha3 (deg)':<12} | {'Alpha4 (deg)':<12} | {'Camber (deg)':<12}")
    print("-" * 62)

    stat_si=stage_index*4+2
    stat_so=stage_index*4+3
    sections=["Hub", "25%", "Mean", "75%", "Tip"]
    
    for i in range(5):
        r_val=rHub+(d['rTip'][stage_index]-rHub)*i/4.0
        a3=np.degrees(d['alpha'][i,stat_si])
        a4=np.degrees(d['alpha'][i,stat_so])
        camber=a3-a4
        print(f"{sections[i]:<10} | {r_val:<8.4f} | {a3:<12.2f} | {a4:<12.2f} | {camber:<12.2f}")


def print_spool(d, pInTot_spool, TInTot_spool):
    stageNum=d['stageNum']
    label=d['label']
    print(f"\n=== {label} MEAN LINE RESULTS ===")
    for s in range(stageNum):
        plane_in=2*s
        plane_rot=2*s+1
        plane_out=2*s+2
        stat_ri=s*4
        stat_ro=s*4+1
        stat_so=s*4+3
        
        V1=d['Vel'][2,stat_ri]
        V2=d['Vel'][2,stat_ro]
        V3=d['Vel'][2,stat_so]
        W1=d['W'][2,stat_ri]
        W2=d['W'][2,stat_ro]
        Vtheta1=d['VelTheta'][2,stat_ri]
        Vtheta2=d['VelTheta'][2,stat_ro]
        
        print(f"\nStage {s+1}")
        print(f"rTip: {d['rTip'][s]:.4f} m, rMean: {d['rMean'][s]:.4f} m")
        print(f"beta1: {np.degrees(d['beta'][2,stat_ri]):.2f} deg, beta2: {np.degrees(d['beta'][2,stat_ro]):.2f} deg")
        print(f"P02: {d['pTot'][2,plane_rot]:.2f} Pa, T02: {d['TTot'][2,plane_rot]:.2f} K")
        
        print(f"--- Velocity Changes ---")
        print(f"W1: {W1:.2f} m/s | W2: {W2:.2f} m/s | Rotor dW (W1-W2): {W1-W2:.2f} m/s")
        print(f"V2: {V2:.2f} m/s | V3: {V3:.2f} m/s | Stator dV (V2-V3): {V2-V3:.2f} m/s")
        print(f"Vtheta1: {Vtheta1:.2f} m/s | Vtheta2: {Vtheta2:.2f} m/s | dVtheta: {Vtheta2-Vtheta1:.2f} m/s")
        
    overall_PR=d['pTot'][2,stageNum*2]/pInTot_spool
    print(f"\nOverall {label} PR: {overall_PR:.4f}")


# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Spool Analysis
lpc=analyse_spool(stageNumLow,omegaLow,pInTot,TempInTot,beta2_mean_LPC,'LPC')

pTot_hpc_in=lpc['pTot'][2,stageNumLow*2]
TTot_hpc_in=lpc['TTot'][2,stageNumLow*2]

hpc=analyse_spool(stageNumHigh,omegaHigh,pTot_hpc_in,TTot_hpc_in,beta2_mean_HPC,'HPC')

# 2. Console Print Outputs
print_spool(lpc,pInTot,TempInTot)
print_spool(hpc,pTot_hpc_in,TTot_hpc_in)

# 3. Overall Results
overall_PR_total=hpc['pTot'][2,stageNumHigh*2]/pInTot
print(f"\n=== OVERALL COMPRESSOR ===")
print(f"Achieved OVP: {overall_PR_total:.2f}")

# 4. CAD Export Prints
print_geometry_for_cad(lpc, 0)
print_geometry_for_cad(hpc, 0)

print_stator_geometry_for_cad(lpc, 0)
print_stator_geometry_for_cad(hpc, 0)
