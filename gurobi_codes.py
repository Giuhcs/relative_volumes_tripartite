# This file is a library containing the functions used to implent optimization tasks with Gurobi along the paper.
# The following packages are needed to run them: gurobipy, numpy.
# Further parametrization following the Gurobi package should be done by hand in the code below.


##### Bilocal scenario #####

# This functions tests whether a given behaviour "p" is local or not by default, returning its non-local distance if not.
# However, by setting test = "bilocality", one can test bilocality, returning its non-bilocal distance if dist = True;
# if this is the case, then gap = 0.1 can be adjusted to set the gap desired between primal and dual solutions to stop the
# optimization. Feasibility (i.e. "dist = False") is always faster than minimization/maximization, so we suggest to set it 
# true if needed only. Besides, by changing tol (10**-6 by default; must be >= 10**-9) one can adjust the feasibility tolerance
# for the primal solution.

def b_compatibility(p,test="locality",dist=False,gap=0.1,tol=10**-6): # p is a list following the indexing p[x*16+z*8+a*4+b*2+c]
	## creating the model ##
	m = gp.Model("local")
	t = m.addVars(range(32),lb=0.0,ub=1.0,vtype=GRB.CONTINUOUS,name="t") # slack variables to compute the distance
	l = m.addVars(range(32),lb=0.0,ub=1.0,vtype=GRB.CONTINUOUS,name="l") # hidden variables indexed as a0*16+a1*8+b*4+c0*2+c1

	## add constraint of hidden distribution summing up to 1 ##
	proper_prob = m.addConstr(sum(l[i] for i in range(32)) == 1,name="pprob")

	## setting the objective ##
	# behaviour derived from hidden distribution l
	q_00 = [sum(l[a0*16+a1*8+b*4+c0*2+c1] for a1 in [0,1] for c1 in [0,1]) for a0 in [0,1] for b in [0,1] for c0 in [0,1]]
	q_01 = [sum(l[a0*16+a1*8+b*4+c0*2+c1] for a1 in [0,1] for c0 in [0,1]) for a0 in [0,1] for b in [0,1] for c1 in [0,1]]
	q_10 = [sum(l[a0*16+a1*8+b*4+c0*2+c1] for a0 in [0,1] for c1 in [0,1]) for a1 in [0,1] for b in [0,1] for c0 in [0,1]]
	q_11 = [sum(l[a0*16+a1*8+b*4+c0*2+c1] for a0 in [0,1] for c0 in [0,1]) for a1 in [0,1] for b in [0,1] for c1 in [0,1]]
	q = q_00+q_01+q_10+q_11

	# constraints from |p-q|<= t
	local_const1 = m.addConstrs(((q[i]-p[i]-t[i]) <= 0 for i in range(32)),name = 'lc1')
	local_const2 = m.addConstrs(((q[i]-p[i]+t[i]) >= 0 for i in range(32)),name = 'lc2')
	
	if test == "bilocality" and dist == False:
		m.setObjective(1,GRB.MINIMIZE)
	else: 
		m.setObjective(sum(t[i] for i in range(32)),GRB.MINIMIZE)
	
	## quadratic constraints for separability between parties A and C ##
	if test == "bilocality":
		q_ac = m.addVars(range(16),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="qac") # auxiliar, following the indexing a0*8+a1*4+c0*2+c1
		q_a = m.addVars(range(4),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="qa") # auxiliar, following the indexing a0*2+a1
		q_c = m.addVars(range(4),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="qc") # auxiliar, following the indexing c0*2+c1

		# correct marginals
		for a0 in [0,1]:
			for a1 in [0,1]:
				for c0 in [0,1]:
					for c1 in [0,1]:
						m.addConstr(q_ac[a0*8+a1*4+c0*2+c1] == sum(q[a0*16+a1*8+b*4+c0*2+c1] for b in [0,1]))
		for a0 in [0,1]:
			for a1 in [0,1]:
				m.addConstr(q_a[a0*2+a1] == sum(q[a0*16+a1*8+b*4+c0*2+c1] for b in [0,1] for c0 in [0,1] for c1 in [0,1]))

		for c0 in [0,1]:
			for c1 in [0,1]:
				m.addConstr(q_c[c0*2+c1] == sum(q[a0*16+a1*8+b*4+c0*2+c1] for b in [0,1] for a0 in [0,1] for a1 in [0,1]))
	
		# separability
		for a0 in [0,1]:
			for a1 in [0,1]:
				for c0 in [0,1]:
					for c1 in [0,1]:
						m.addQConstr(q_ac[a0*8+a1*4+c0*2+c1] == q_a[a0*2+a1]*q_c[c0*2+c1])

	## setting optimization parameters, updating and solving the model ##
	m.setParam('OutputFlag', False) # this is to make the output quiet
	if test == "bilocality":
		m.params.NonConvex = 2 # this tells the model to consider a non-convex feasibility region; with quadratic (2) constraints
		if dist == True:
			m.setParam('MIPGap',gap)
	m.setParam('FeasibilityTol',tol) # set tolerance for assessing primal feasibility
	m.update()
	m.optimize()

	out = 0
	if test == "bilocality" and dist == True:
		out = m.ObjVal
	else:
		out = m.Status
	return out

##### Evans scenario #####

# This function tests whether a given behaviour "beh" is compatible or not with the classical Evans scenario,
# returning its non-evans distance if #dist = True"; if this is the case, then "gap = 0.1" can be adjusted to set the gap desired between 
# primal and dual solutions to stop the optimization. Feasibility (i.e. "dist = False") is always faster than 
# minimization/maximization, so we suggest to set it true if needed only. Besides, by changing tol (10**-6 by default; must be >= 10**-9)
# one can adjust the feasibility tolerance for the primal solution. It is possible to test compatibility with a generic p(a,b,c) by setting
# from_simplex = True. Finally, by setting intervention = True one can also take into account constraints from intervention.

def e_compatibility(beh, from_simplex=False, intervention=False,dist=False, gap=0.1, tol=10**-6): # beh is a list following the indexing beh[x*16+z*8+a*4+b*2+c]
	## creating the model ##
	m = gp.Model()
	q = m.addVars(range(32),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="q") # hidden variables, following the indexing a0*16+a1*8+b*4+c0*2+c1
	if dist == True:
		t = m.addVars(range(8),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="t") # slack variables to compute the distance
	
	## constraint of being a proper joint probability distribution ##
	m.addConstr(sum(q[i] for i in range(32)) == 1)
	
	## quadratic constraints for separability between parties A and C ##
	q_ac = m.addVars(range(16),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="qac") # auxiliar, following the indexing a0*8+a1*4+c0*2+c1
	q_a = m.addVars(range(4),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="qa") # auxiliar, following the indexing a0*2+a1
	q_c = m.addVars(range(4),lb=0,ub=1,vtype=GRB.CONTINUOUS,name="qc") # auxiliar, following the indexing c0*2+c1

	# correct marginals
	for a0 in [0,1]:
		for a1 in [0,1]:
			for c0 in [0,1]:
				for c1 in [0,1]:
					m.addConstr(q_ac[a0*8+a1*4+c0*2+c1] == sum(q[a0*16+a1*8+b*4+c0*2+c1] for b in [0,1]))
	for a0 in [0,1]:
		for a1 in [0,1]:
			m.addConstr(q_a[a0*2+a1] == sum(q[a0*16+a1*8+b*4+c0*2+c1] for b in [0,1] for c0 in [0,1] for c1 in [0,1]))

	for c0 in [0,1]:
		for c1 in [0,1]:
			m.addConstr(q_c[c0*2+c1] == sum(q[a0*16+a1*8+b*4+c0*2+c1] for b in [0,1] for a0 in [0,1] for a1 in [0,1]))
	
	# separability
	for a0 in [0,1]:
		for a1 in [0,1]:
			for c0 in [0,1]:
				for c1 in [0,1]:
					m.addQConstr(q_ac[a0*8+a1*4+c0*2+c1] == q_a[a0*2+a1]*q_c[c0*2+c1])

	## setting the objective ##
	if from_simplex = False:
		p = [beh[b*16+b*8+a*4+b*2+c] for a in [0,1] for b in [0,1] for c in [0,1]] # projected behaviour
	else:
		p = beh
	p_q = [sum(q[a*16+a1*8+b*4+c*2+c1] for a1 in [0,1] for c1 in [0,1]) if b == 0 else sum(q[a0*16+a*8+b*4+c0*2+c] for a0 in [0,1] for c0 in [0,1]) for a in[0,1] for b in [0,1] for c in [0,1]] # behaviour p_q from hidden distribution q
	if dist == False:
		# constraint for producing the correct projected behaviour
		for a in [0,1]:
			for b in [0,1]:
				for c in [0,1]:
					m.addConstr(p[a*4+b*2+c] == p_q[a*4+b*2+c])
		m.setObjective(1,GRB.MINIMIZE)
	else:
		# constraints from |p-p_q|<= t
		for i in range(8):
			m.addConstr(p[i]-p_q[i]-t[i] <= 0)
			m.addConstr(p[i]-p_q[i]+t[i] >= 0)
		m.setObjective(sum(t[i] for i in range(8)),GRB.MINIMIZE)

	## constraints from intervention ##
	# constraints from do b=0
	for a in [0,1]:
		m.addConstr(sum(q_a[a*2+a1] for a1 in [0,1]) == sum(p[a*4+b*2+c] for b in [0,1] for c in [0,1]))
	for c in [0,1]:
		m.addConstr(sum(q_c[c*2+c1] for c1 in [0,1]) == sum(p[a*4+b*2+c] for a in [0,1] for b in [0,1]))

	# constraints from do b=1
	for a in [0,1]:
		m.addConstr(sum(q_a[a0*2+a] for a0 in [0,1]) == sum(p[a*4+b*2+c] for b in [0,1] for c in [0,1]))
	for c in [0,1]:
		m.addConstr(sum(q_c[c0*2+c] for c0 in [0,1]) == sum(p[a*4+b*2+c] for a in [0,1] for b in [0,1]))

	## setting optimization parameters, updating and solving the model ##
	m.setParam('OutputFlag', False) # this is to make the output quiet
	m.params.NonConvex = 2 # this tells the model to consider non-convex feasibility region; with quadratic (2) constraints
	m.setParam('FeasibilityTol',10**-6) # set tolerance for assessing primal feasibility
	m.setParam('MIPGap',gap)
	m.update()
	m.optimize()
	
	out = 0
	if dist == False:
		out = m.Status
	else:
		out = m.ObjVal
	return out

# This function tests whether a given behaviour "p" is NSI, be checking if there is a NSI behaviour from the bilocal scenario
# that generates it

def e_NSI_test(p):
	m=gp.Model()
	P=m.addMVar((2,2,2,2,2),lb=0,ub=1)
	#Normalization of P(a,b,c|x,z)
	m.addConstrs(sum(P[a,b,c,x,z]for a in range(2)for b in range(2)for c in range(2))==1 for x in range(2)for z in range(2))

	#Non-signaling P(a|x,z)=P(a|x), P(c|x,z)=P(c|z), P(b|x,z)=P(b)
	m.addConstrs(sum(P[a,b,c,x,0]for b in range(2)for c in range(2))==sum(P[a,b,c,x,1]for b in range(2)for c in range(2))for a in range(2)for x in range(2))
	m.addConstrs(sum(P[a,b,c,0,z]for b in range(2)for a in range(2))==sum(P[a,b,c,1,z]for b in range(2)for a in range(2))for c in range(2)for z in range(2))
	m.addConstrs(sum(P[a,b,c,0,0]for a in range(2)for c in range(2))==sum(P[a,b,c,x,z]for a in range(2)for c in range(2))for b in range(2)for x in range(2)for z in range(2))

	#d-separation p(a,c|x,z)=p(a|x)p(c|z)
	m.addConstrs(sum(P[a,b,c,x,z]for b in range(2))==sum(P[a,b1,c_,x,0]*P[a_,b2,c,0,z] for b1 in range(2) for c_ in range(2)for a_ in range(2) for b2 in range(2)) for a in range(2) for c in range(2) for x in range(2) for z in range(2))

	#Compatibility constraints
	m.addConstrs(P[a,b,c,b,b]==p[a,b,c] for a in range(2) for b in range(2) for c in range(2))

	m.Params.OutputFlag=0
	m.Params.MIPGap=1e-6
	m.Params.OptimalityTol=1e-6
	m.params.NonConvex = 2
	m.params.TimeLimit = 1

	m.setObjective(1, gp.GRB.MINIMIZE)
	m.optimize()

	return m.Status