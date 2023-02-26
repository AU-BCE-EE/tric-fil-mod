
def Kga_onda(pH, TK, KH, pKa, P, ssa, Qstd, Qliq, por, Ax): 
    
    dp = 6*(1-por)/ssa  # characteristic packing length

    if dp < 15:
        dp_emp = 2.0
    else:
        dp_emp = 5.23

    Q = Qstd/3600
    v_g = Q/Ax  # linear velocity
    v_l = 1.59e-5   # m/sec
    rho_g = P*28.97/(0.08206*TK)*1000
    visc_g = 9.1e-8*TK-1.16e-5   # empirical relation for gas viscosity vs TK
    Dg = 1.16e-5    # gas diffusion coefficient in m2/sec; compound specific
    visc_l = -2.55e-5*TK+8.51e-3
    rho_l = 1.0e6   # g/m3
    Dliq = 1.89e-9  # liquid diffusion coefficient
    g = 9.81    # m/sec^2
    sigm_c = 0.75   # critical surface tension
    sigm_l = 0.0073 # surface tension
    Re = rho_l*v_l/(ssa*visc_l)
    Fr = v_l*v_l*ssa/g
    We = v_l*v_l*rho_l/(sigm_l*ssa)
      
    ae = ssa*(1.0-2.71828**(-1.45*(sigm_c/sigm_l)**0.75 *
                            Re**0.1*Fr**-0.05*We**0.2))
    
    kg = dp_emp*(v_g*rho_g/(ssa*visc_g))**0.7 \
        * (visc_g/(rho_g*Dg))**(1/3) * (ssa*dp)**-2 * ssa*Dg
    
    kl = 0.0051*(v_l*rho_l/(
        ae * visc_l))**(2/3) * (visc_l/(rho_l * Dliq))**(-0.5) * (
        ssa*dp)**0.4 * (rho_l/(visc_l*g))**(-1/3)
    
    R = 0.083144  # Gas constant (L bar / K-mol)
    Kaw = KH/(R*TK) #  Neutral air-water distribution
    alpha = 1/(1+10**(pH - pKa))
    Daw = alpha*Kaw
    
    Rtot = 1/(kg*ae) + Daw/(kl*ae)
    Kga = 1/Rtot
    return Kga
