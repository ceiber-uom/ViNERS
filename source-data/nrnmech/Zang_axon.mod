TITLE Zang Axon Internode channels

: 09/21
: Calvin Eiber
:
: Merged from .MOD files nadp, leak
:
: 08/21
: Yunliang Zang (as uploaded to modelDB)
:
: Includes the following: 
:  leak  the leak current
:  nadk  Sodium ion accumulation with radial and longitudinal diffusion and pumping
: 		 initial Na is 10 mM NEURON provided default, for the present strategy, the 
: 		 value of k4 makes initial contration equilibrates at this value.
:  

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)

    (mol) = (1)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    FARADAY = (faraday) (10000 coulomb)
    R = (k-mole) (joule/degC)    
    PI = (pi) (1)
}

NEURON {
    SUFFIX zang_axon
    NONSPECIFIC_CURRENT il
    USEION na READ nai, ina WRITE nai, ina, nao
    
    NONSPECIFIC_CURRENT ik_pump
    
    RANGE gbar_leak
    RANGE ina_pump, TotalPump, ik_ratio, na, pump, pumpna, na3,DNa, k1, k2, k3, k4, pcu,nai_i
    GLOBAL vrat
}

DEFINE Nannuli 4

CONSTANT { volo = 1e14 (um2) }

PARAMETER {
    gbar_leak  = 0.000125 (S/cm2)
	Eleak = -65 (mV)

	nai_i = 10
	nao_i = 140
    DNa = 0.6 (um2/ms)
    k1 = 2.0 (/mM3-ms)
    k2 = 0.001 (/ms)
    k3 = 0.6 (/ms)
                            : to eliminate pump, set TotalPump to 0 in hoc
    TotalPump = 1e-14 (mol/cm2)
    
    ik_ratio = -0.66666666 (1)
}

ASSIGNED { 
    v (mV)
    il (mA/cm2)

    diam (um)
    L (um)
    ina (mA/cm2)
    nai (mM)
    vrat[Nannuli] : numeric value of vrat[i] equals the volume
                  : of annulus i of a 1um diameter cylinder
                  : multiply by diam^2 to get volume per um length
    
    k4          (/mM3-ms)

    ina_pump (mA/cm2)
    parea (um)

    ik_pump (mA/cm2)
	pcu (mA/cm2)	
}

STATE {
    na[Nannuli] (mM) <1e-3>
    nao (mM)
    pump (mol/cm2)
    pumpna (mol/cm2)
}

BREAKPOINT {

    SOLVE state METHOD sparse

    ina = ina_pump
    ik_pump = ik_ratio*ina_pump
    pcu = ik_pump + ina_pump

    il = gbar_leak*(v-Eleak)
}

LOCAL factors_done

INITIAL {
	nai = nai_i
	nao = nao_i
    k4=(((10/nao)^3)*k1*k3)/k2    :Set the equilibrium at nai0_na_ion
    parea = PI*diam
    pump = TotalPump/(1 + (nai*k1/k2))
    pumpna = TotalPump - pump
    if (factors_done == 0) {    : flag becomes 1 in the first segment
        factors_done = 1        : all subsequent segments will have
        factors()               : vrat = 0 unless vrat is GLOBAL
    }

    FROM i=0 TO Nannuli-1 {
        na[i] = nai
    }
}

LOCAL frat[Nannuli]     : scales the rate constants for model geometry

PROCEDURE factors() {
    LOCAL r, dr2
    r = 1/2                 : starts at edge (half diam)
    dr2 = r/(Nannuli-1)/2   : full thickness of outermost annulus,
                            : half thickness of all other annuli
    vrat[0] = 0
    frat[0] = 2*r
    FROM i=0 TO Nannuli-2 {
        vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  : interior half
        r = r - dr2
        frat[i+1] = 2*PI*r/(2*dr2)              : outer radius of annulus
                                                : div by distance between centers
        r = r - dr2
        vrat[i+1] = PI*(r+dr2/2)*2*dr2 : outer half of annulus
    }
}

KINETIC state {
    COMPARTMENT i, diam*diam*vrat[i] {na}
    COMPARTMENT (1e10)*parea {pump pumpna}
    COMPARTMENT volo {nao}

    LONGITUDINAL_DIFFUSION i, DNa*diam*diam*vrat[i] {na}

    :pump
    ~ 3 na[0] + pump <-> pumpna (k1*parea*(1e10), k2*parea*(1e10))
    ~ pumpna <-> pump + 3 nao (k3*parea*(1e10), k4*parea*(1e10))
    ina_pump = FARADAY*(f_flux - b_flux)/parea
    CONSERVE pump + pumpna = TotalPump * parea * (1e10)

    : all currents except pump
    ~ na[0] << (-(ina-ina_pump)*PI*diam/(FARADAY)) : ina is Na efflux, subtracting ina_pump this is quite right!
	
    FROM i=0 TO Nannuli-2 {
        ~ na[i] <-> na[i+1] (DNa*frat[i+1], DNa*frat[i+1])
    }

    nai = na[0]
}
