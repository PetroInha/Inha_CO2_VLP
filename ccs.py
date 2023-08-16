import numpy as np
from math import *

def density_of_co2(P, T, verbose = False):
    """ this is to calculate density of co2

    Args:
        P (float): Pressure [bar]
        T (float): Temperature [K]
        verbose (bool, optional): Whether to print out result. Defaults to False.

    Returns:
        density (float): calculated density of CO2 [kg/m^3]
    """
    if P>25 and P<=100:
        A1 = 2.089800972761597e5;   B1 = -1.456286332143609e4; C1 = 2.885813588280259e2;   D1 = -1.597103845187521
        A2 = -1.675182353338921e3;  B2 = 1.16799554255704e2;   C2 = -2.31558333122805;     D2 = 1.284012022012305e-2
        A3 = 4.450600950630782;     B3 = -3.10430147581379e-1; C3 = 6.157718845508209e-3;  D3 = -3.420339567335051e-5
        A4 = -3.919844561756813e-3; B4 = 2.734973744483903e-4; C4 = -5.428007373890436e-6; D4 = 3.019572090945029e-8
        
        alpha = A1 + B1*P + C1*P**2 + D1*P**3
        beta  = A2 + B2*P + C2*P**2 + D2*P**3
        gamma = A3 + B3*P + C3*P**2 + D3*P**3
        theta = A4 + B4*P + C4*P**2 + D4*P**3

        density = alpha + beta*T + gamma*(T**2) + theta*(T**3)
        if verbose == True:
            print('density: {0}'.format(density))
    
    elif P>100 and P<700:
        A1 = 1.053293651041897e5;   B1 = -9.396448507019846e2;  C1 = 2.397414334181339;     D1 = -1.819046028481314e-3
        A2 = -8.253383504614545e2;  B2 = 7.618125848567747;     C2 = -1.963563757655062e-2; D2 = 1.497658394413360e-5
        A3 = 2.135712083402950;     B3 = -2.023128850373911e-2; C3 = 5.272125417813041e-5;  D3 = -4.043564072108339e-8
        A4 = -1.827956524285481e-3; B4 = 1.768297712855951e-5;  C4 = -4.653377143658811e-8; D4 = 3.586708189749551e-11

        alpha = A1 + B1*P + C1*P**2 + D1*P**3
        beta  = A2 + B2*P + C2*P**2 + D2*P**3
        gamma = A3 + B3*P + C3*P**2 + D3*P**3
        theta = A4 + B4*P + C4*P**2 + D4*P**3

        density = alpha + beta*T + gamma*(T**2) + theta*(T**3)
        if verbose == True:
            print('density: {0}'.format(density))

    else:
        assert False, "Warring: co2 pressure is not valid"
    
    return density

def viscosity_of_co2(density, T, verbose = False):
    """this is to calculate viscosity of co2

    Args:
        density (float): Density [kg/m^3]
        T (float): Temperature [K]
        verbose (bool, optional): Whether to print out result. Defaults to False.

    Returns:
        viscosity (float): calculated viscosity of CO2 [Pa s]
    """
    a0 = 0.235156; a1 = -0.491266; a2 = 5.211155e-2; a3 = 5.347906e-2; a4 = -1.537102e-2
    d11 = 0.4071119e-2; d21 = 0.7198037e-4; d64 = 0.2411697e-16; d81 = 0.2971072e-22 ; d82 = -0.1627888e-22
    
    viscosity1 = 1.00697*sqrt(T)/np.exp(a0 + a1*log(T/251.196) + a2*(log(T/251.196))**2 + a3*(log(T/251.196))**3 + a4*(log(T/251.196))**4)   
    viscosity2 = d11*density + d21*(density**2) + d64*(density**6)/((T/251.196)**3) + d81*(density**8) + d82*(density**8)/(T/251.196)

    viscosity = viscosity1 + viscosity2

    if verbose == True: 
        print('viscosity: {0}'.format(viscosity))

    return viscosity*1E-6

def Reynolds(density, V, D, viscosity):
    """this is to calculate Reynolds number of co2

    Args:
        density (float): Density [kg/m^3]
        V (float): Velocity [m/s]
        D (float): Diameter [m]
        viscosity (float): [Pa s]

    Returns:
        Reynolds (float): calculated Reynolds number of co2 [dimensionless]
    """
    Re = density*V*D/(viscosity)
    return Re

def friction_factor(rel_roughness, Re):
    """this is to calculate friction factor of CO2

    Args:
        rel_roughness (float): relative roughness [dimensionless]
        Re (float): Reynolds number [dimension less]

    Returns:
        friction_factor(float) : calculated friction factor [dimensionless]
    """
    if Re > 4000:
        AA0 = -0.79638*log(rel_roughness/8.208 + 7.3357/Re)
        AA1 = Re*rel_roughness + 9.3120665*AA0
        f = ((8.128943 + AA1)/(8.128943*AA0) + 0.86859509*AA1*log(AA1/(3.7099535*Re)))
    elif Re < 2100:
        f = 64/Re
    else:
        f = np.nan
    return f

def flow_rate(M,density):
    """this is to calculate volume flow rate of CO2

    Args:
        M (float): mass flow rate [kg/s]
        density (float): density [kg/m^3]

    Returns:
        flow_rate(float) : calculated volume flow rate [m^3/s]
    """
    Q=M/density
    return Q

def velocity(D,Q):
    """this is to calculate velocity of CO2

    Args:
        D (float): diameter [m]
        Q (float): flow rate [m^3/s]

    Returns:
        velocity(float): calculated velocity [m/s]
    """
    V = Q/(pi*(D**2)/4)
    return V

def friction_delta_pressure(f,viscosity,V,D,h):
    """this is to calculate pressure loss due to friction

    Args:
        f (float): friction factor [dimensionless]
        viscosity (float): viscosity [Pa s]
        V (float): velocity [m/s]
        D (float): diameter [m]
        h (float): height of pipe [m]

    Returns:
        friction_delta_pressure: calculated friction delta pressure [Pa]
    """
    delta_P = f*viscosity*h*(V**2)/(2*D)
    return delta_P

def compute_next_pressure(P1, T, D, M, length, del_depth,rel_roughness, gravity_acc = 9.81, bar2pa = 1E5, Accepted_error = 1E-6):
    """ this is to compute next node's pressure

    Args:
        P1 (float): current node pressure
        T (float): temperature
        D (float): diameter
        M (float): CO2 mass flow rate 
        length (float): flow length between two nodes
        del_depth (float): diff of depth between two nodes

    Returns:
        P2 (float): next node pressure
    """
    # Guess only by gravity 
    assumed_co2_density = density_of_co2(P1, T)
    hydro_static_p = assumed_co2_density*gravity_acc*del_depth/bar2pa
    P2 = P1 + hydro_static_p
    P_avg = (P1 + P2)/ 2

    # Iteration to compute next pressure 
    error = 99999
    while error > Accepted_error:
        
        # CO2 - EOS 
        assumed_co2_density = density_of_co2(P_avg, T)
        assumed_co2_viscosity = viscosity_of_co2(assumed_co2_density, T)
        
        # To get pressure change due to friction
        q = flow_rate(M, assumed_co2_density)
        flow_velocity = velocity(D, q)
        reynolds = Reynolds(assumed_co2_density, flow_velocity, D, assumed_co2_viscosity )
        f = friction_factor(rel_roughness, reynolds )
        P_del_friction = friction_delta_pressure(f,assumed_co2_viscosity,flow_velocity,D,length)
        
        # To get pressrue change due to gravity 
        P_del_gravity = assumed_co2_density*gravity_acc*del_depth/bar2pa

        P2_new = P1 + P_del_gravity - P_del_friction/bar2pa

        error = P2 - P2_new
        # print(f'error: {error}')

        P2 = P2_new
        P_avg = (P1 + P2)/ 2
        
    return P2

if __name__ == "__main__":


    P_up = input('압력(bar): ')
    T = 300.0              # k
    h = 100.0              # m
    D = 0.25               # m
    M = 2.87666            # kg/s
    rel_roughness = 0.0018

    P = (P_up + P_up+2.0)/2.

    density = density_of_co2(P, T)
    viscosity_of_co2 = viscosity_of_co2(density, T)

    viscosity = viscosity_of_co2

    Q = flow_rate(M,density)

    V = velocity(float(D),float(Q))

    Re = Reynolds(density, float(V), float(D), viscosity)

    rel_roughness = float(rel_roughness)

    friction_factor(float(rel_roughness),Re)
    f = friction_factor(float(rel_roughness),Re)

    def pressure_under(P,delta_P):
        P_under = P + delta_P
        print('입력한 압력보다 100m아래의 압력은 {0}(bar)입니다.'.format(P_under))
