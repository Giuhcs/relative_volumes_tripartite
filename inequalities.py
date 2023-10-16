# This is a library containing the functions to test compatibility of a given p(a,b,c) with the triangle scenario by checking inequalities.
# Required packages: numpy, math

#### Finner inequality ####

#def finner_test(p):
#	return 0

#### Entropic inequalities ####

# entropy h(p) of a distribution p
def h(p): # here p is a list of 8 entries
	aux = 0
	for i in range(len(p)):
		aux = aux - p[i]*math.log(p[i])
	return aux

# entropic inequalities
def ineq_1(p_abc,p_a,p_b,p_c):
	return 2*h(p_abc) >= h(p_a) + h(p_b) + h(p_c)

def ineq_21(p_ab,p_bc,p_a,p_b,p_c):
	return h(p_ab) + h(p_bc) >= h(p_a) + h(p_b) + h(p_c)

def ineq_22(p_bc,p_ac,p_a,p_b,p_c):
	return h(p_bc) + h(p_ac) >= h(p_a) + h(p_b) + h(p_c)

def ineq_23(p_ac,p_ab,p_a,p_b,p_c):
	return h(p_ac) + h(p_ab) >= h(p_a) + h(p_b) + h(p_c)

def ineq_31(p_ab,p_c,p_bc,p_ac):
	return h(p_ab) + h(p_c) <= h(p_bc) + h(p_ac)

def ineq_32(p_bc,p_a,p_ab,p_ac):
	return h(p_bc) + h(p_a) <= h(p_ab) + h(p_ac)

def ineq_33(p_ac,p_b,p_bc,p_ab):
	return h(p_ac) + h(p_b) <= h(p_bc) + h(p_ab)

# main function; return 0 if p satisfies all the inequalities listed above and 1 otherwise
def entropic_test(p):
	p_ab = [sum(p[a*4+b*2+c] for c in [0,1]) for a in [0,1] for b in [0,1]]
	p_ac = [sum(p[a*4+b*2+c] for b in [0,1]) for a in [0,1] for c in [0,1]]
	p_bc = [sum(p[a*4+b*2+c] for a in [0,1]) for b in [0,1] for c in [0,1]]
	p_a = [sum(p[a*4+b*2+c] for b in [0,1] for c in [0,1]) for a in [0,1]]
	p_b = [sum(p[a*4+b*2+c] for a in [0,1] for c in [0,1]) for b in [0,1]]
	p_c = [sum(p[a*4+b*2+c] for a in [0,1] for b in [0,1]) for c in [0,1]]
	out = 1
	if ineq_1(p,p_a,p_b,p_c) & ineq_21(p_ab,p_bc,p_a,p_b,p_c) & ineq_22(p_bc,p_ac,p_a,p_b,p_c) & ineq_23(p_ac,p_ab,p_a,p_b,p_c) & ineq_31(p_ab,p_c,p_bc,p_ac) & ineq_32(p_bc,p_a,p_ab,p_ac) & ineq_33(p_ac,p_b,p_bc,p_ab):
		out = 0
	return out