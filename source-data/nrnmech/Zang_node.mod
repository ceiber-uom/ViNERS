TITLE Zang Axon Node channels

: 09/21
: Calvin Eiber
:
: Merged from .MOD files ka, kd, na
:
: 08/21
: Yunliang Zang (as uploaded to modelDB)
:
: Includes the following: 
:  NA  the transient, inward sodium current
:  KA  the transient, outward potassium current   // change it to Connor-Stevens model
:  KD  the delayed, rectifying potassium current  // change it to Connor-Stevens model
:
:  the other node channels are shared, see zang_internode.mod


NEURON {
	SUFFIX zang_node

    
    USEION na READ ena WRITE ina
    USEION k READ ko WRITE ik

	RANGE gbar_na, ina, q10na   : SUFFIX na
    RANGE gbar_ka, ika, q10ka   : SUFFIX ka    
    RANGE gbar_kd, ikd, q10kd	: SUFFIX kd

    THREADSAFE
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}


CONSTANT {
	FARADAY = 96485.3399
	R = 8.31447215
}



PARAMETER {

	: NA
	gbar_na = 0.014 (S/cm2)
	q10na = 1	

	: KA
	gbar_ka  = 0.005 (S/cm2)	
    q10ka = 1

	: KD
    gbar_kd = 0.003 (S/cm2)	
	q10kd = 1

	Ekd1 = -70 (mV)
}

ASSIGNED {

	: NA
    v (mV)
    ina (mA/cm2)
    ena
	
	mna_inf
	tau_mna
	hna_inf
	tau_hna
	qtna
	
	:KA
	ko
	ik  (mA/cm2)
	ika (mA/cm2)
	ikd (mA/cm2)

	mka_inf
	tau_mka
	hka_inf
	tau_hka
	qtka

	:KD
	am
	bm
	mkd_inf
	tau_mkd
	qtkd

    celsius (degC)
}

STATE {
 	mna hna mka hka mkd
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gbar_na*mna*mna*mna*hna*(v-ena)  : NA
    ika = gbar_ka*mka*mka*mka*hka*(v-Ekd1) : KA
	ikd = gbar_kd*mkd*mkd*mkd*mkd*(v-Ekd1) : KD
	ik = ika+ikd
}

INITIAL {
	settables(v)
	
    mna = mna_inf
    hna = hna_inf
    qtna = q10na^((celsius-6.3 (degC))/10 (degC))

    mka = mka_inf
    hka = hka_inf
    qtka = q10ka^((celsius-6.3 (degC))/10 (degC))

    mkd = mkd_inf
    qtkd = q10kd^((celsius-6.3 (degC))/10 (degC))
}

DERIVATIVE states {
	settables(v)
	
    mna' = (mna_inf-mna)/tau_mna
    hna' = (hna_inf-hna)/tau_hna

    mka' = (mka_inf-mka)/tau_mka
    hka' = (hka_inf-hka)/tau_hka

    mkd' = (mkd_inf-mkd)/tau_mkd
}

UNITSOFF

PROCEDURE settables(v (mV)) { 
	
	: NA
	mna_inf = 1/(1+exp(-(v+48-10)/8.5))
	tau_mna = (0.132/cosh((v+27)/7.5)+0.003/(1+exp(-(v+27)/5)))/qtna
	hna_inf = 1/(1+exp((v+47)/6))
	tau_hna = 10/cosh((v+42)/15)/qtna

	: KA
	mka_inf = (.0761*exp((v+94.22)/31.84)/(1+exp((v+1.17)/28.93)))^(.3333)
	tau_mka = (.3632+1.158/(1+exp((v+55.96)/20.12)))/qtka
	hka_inf =1/(1+exp((v+53.3)/14.54))^4
	tau_hka = (1.24+2.678/(1+exp((v+50)/16.027)))/qtka

	: KD
	am = (-.01*(v+50-4.3)/(exp(-(v+50-4.3)/10)-1))/2*3.8
	bm = .125*exp(-(v+60-4.3)/80)/2*3.8
	mkd_inf = am/(am+bm)
	tau_mkd = 1/(am+bm)*3.8/qtkd
}

UNITSON