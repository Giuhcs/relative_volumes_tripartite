# This file is a library containing the function used to implent optimization tasks with the Inflation package along the paper.
# The following packages are needed to run it: inflation, numpy.

##### General #####

# This function tests whether a given behaviour "p" is compatible or not with the provided "scenario"
# for a given inflation level ("inf_lvl") - meaning the number of instances for each source - and a given relaxation level
# ("npa_lvl") for the NPA routine underlying the test, returning the feasibility status provided by the inflation package. It 
# tests compatibility with classical inflations by default, but compatibility with a quantum inflations are possible by
# setting commuting = False.

def inflation_test(p,scenario,inf_lvl,npa_lvl,commuting=True):
# p is a numpy.ndarray: p[a,b,c,x,y,z], with a,b,c being the cardinalities of the (here binary) outputs for parties 
# A, B and C, respectively, and x, y and z those for the inputs of the same respective parties depending on the scenario;
# scenario is a string among: "bilocal", "evans" and "triangle";
# lvl and npa are integers >= 1.
	if scenario == "bilocal":
		dag = {"rho_AB":["A","B"],"rho_BC":["B","C"]}
		nr_inputs = [2,1,2]
		nr_copies = [inf_lvl,inf_lvl]
	elif scenario == "evans":
		dag = {"rho_AB":["A","B"],"rho_BC":["B","C"],"B":["A","C"]}
		nr_inputs = [1,1,1]
		nr_copies = [inf_lvl,inf_lvl]
	else:
		dag = {"rho_AB":["A","B"],"rho_BC":["B","C"],"rho_AC":["A","C"]}
		nr_inputs = [1,1,1]
		nr_copies = [inf_lvl,inf_lvl,inf_lvl]
	nr_outputs = [2,2,2]
	scenario = InflationProblem(dag,nr_outputs,nr_inputs,nr_copies)
	sdp = InflationSDP(scenario,commuting=commuting)
	sdp.generate_relaxation( "npa{}".format(npa_lvl))
	sdp.set_distribution(p)
	sdp.solve()
	return sdp.status