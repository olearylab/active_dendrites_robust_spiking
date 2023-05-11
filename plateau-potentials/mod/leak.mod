
COMMENT
A constant leak current, intended to simulate long-lived inhibition.
ENDCOMMENT

NEURON {
    SUFFIX leak
    NONSPECIFIC_CURRENT i
    RANGE g, i
    RANGE e
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	g = 0.008 (mS/cm2)
	e = -80	(mV)
}

ASSIGNED {
	v	(mV)
	i	(mA/cm2)
}

BREAKPOINT {
	i = g*(v - e)
}
