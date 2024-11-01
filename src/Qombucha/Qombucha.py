
import gurobipy as gp, numpy as np, math
from numpy import linalg as LA
import sys, os, argparse
from copy import deepcopy


def LP_get_P(C,B,c,t,weight_per_patient):
    #C is the cell type matrix C_i is for cel type i 
  
    n,m = C.shape # n cell type, m features
    r,_m = B.shape # r patients, m features
    
    # assert B and C have the same number of columns
    assert m ==_m
    
    model = gp.Model()
    model.Params.Threads = c
    model.Params.TimeLimit = t

    # initialize variables for inferred fractions for all patients
    P = np.empty((r,n),dtype=object)
    for i in range(r):
        for j in range(n):
            P[i,j]= model.addVar(lb=0.0, ub=1.0, vtype= gp.GRB.CONTINUOUS,name='var_p_{}_{}'.format(i,j))
 
    # set constraints for inferred fractions

    for i in range(r):
        # for j in range(n):
            # model.addConstr(P[i,j]>=0.0)
            # model.addConstr(P[i,j]<=1.0)
        model.addConstr(gp.quicksum(P[i]) ==1.0,name = 'constr_p_{}'.format(i))
   

    # set objective

    #initialize B_hat and store all the entries in it
       # set objective

    B_hat = np.matmul(P,C)
    #residue = np.abs(B- B_hat).flatten()
    residue = B- B_hat

    residue = residue * weight_per_patient[:, np.newaxis]
    
    residue = residue.flatten()

    Obj = gp.quicksum(r**2 for r in residue)
    model.setObjective(Obj, gp.GRB.MINIMIZE)

    model.optimize()

    P_arr = np.empty((r,n),dtype=float)
    # _C_arr = np.empty((n,m),dtype=float)
    
    for i in range(r):
        for j in range(n):
            P_arr[i,j] = P[i,j].X

    return P_arr, model.getObjective().getValue()


def LP_get_C(P,B,_C, c,t, weight_per_patient, Alpha, C_initial):
    
    #C is the cell type matrix C_i is for cel type i 
  
    r_,n = P.shape # n cell type, m features
    r,m = B.shape # r patients, m features
    
    # assert B and P have the same number of rows
    assert r == r_
    
    model = gp.Model()
    model.Params.Threads = c
    model.Params.TimeLimit = t

    # initialize variables for unknown entries in _C
    C = _C.copy().astype(object)
        # initialize variables for inferred entries
    
    for i in range(n):
        for j in range(m): 
            if np.isnan(_C[i,j]):
                C[i,j] = model.addVar(lb=0.0, ub=1.0, vtype= gp.GRB.CONTINUOUS,name='var_c_{}_{}'.format(i,j))

    # set objective

    #initialize B_hat and store all the entries in it
    # set objective

    B_hat = np.matmul(P,C)
    #residue = np.abs(B- B_hat).flatten()
    residue = B- B_hat

    residue = residue * weight_per_patient[:, np.newaxis]
    
    residue = residue.flatten()
    
    Obj = gp.quicksum(r**2 for r in residue)

    if Alpha:
        quad_expr = gp.QuadExpr(0)
        for i in range(n):
            for j in range(m):
                if np.isnan(_C[i,j]):
                    quad_expr.add((C[i,j]-C_initial[i,j])**2, Alpha)
        Obj = gp.quicksum((Obj, quad_expr))

    model.setObjective(Obj, gp.GRB.MINIMIZE)

    model.optimize()

    C_arr = np.empty((n,m),dtype=float)
    
    for i in range(n):
        for j in range(m):
            if np.isnan(_C[i,j]):
                C_arr[i,j] = C[i,j].X
            else:
                C_arr[i,j] = _C[i,j]
    
    return C_arr, model.getObjective().getValue()





def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-C', '--reference', 
        type=str, 
        required=True,
        help='Complete reference panel')
    parser.add_argument(
        '-Cp', '--reference_partial', 
        type=str, 
        required=True,
        help='Partial reference panel')
    parser.add_argument(
        '-B', '--bulk', 
        type=str, 
        required=True,
        help='Bulk methylation value matrix')
    parser.add_argument(
        '-op', '--output_P', 
        type=str, 
        required=True,
        help='Path to output P matrix (cell type fractions per patient)')
    parser.add_argument(
        '-oc', '--output_C', 
        type=str, 
        required=True,
        help='Path to output C matrix (cell types by methylation features)')
    parser.add_argument(
        '-I', '--iters', 
        type=int, 
        required=False, 
        default=200,
        help='Max number of iterations (default = 200)')
    parser.add_argument(
        '-c', '--threads', 
        type=int, 
        required=False, 
        default=2,
        help='Number of threads to use (default = 2)')
    parser.add_argument(
        '-t', '--max_run_time_per_gurobi_instance', 
        type=int, 
        required=False, 
        default=100000,
        help='Max runtime per Gurobi instance in seconds (default = 100000)')
    parser.add_argument(
        '-Delta', '--delta_percent_obj_error_threshold', 
        type=float, 
        required=False, 
        default=0.0125,
        help='Convergence threshold value for relative objective error (default = 0.0125)')
    parser.add_argument(
        '-Lambda', '--lambda_diff_successive_delta', 
        type=float, 
        required=False, 
        default=1e-9,
        help='Convergence threshold value for change in successive deltas (default = 1e-9)')
    parser.add_argument(
        '-W', '--weight_per_patient', 
        type=str, 
        required=False, 
        help='Experimental: Weight per patient to be used in the objective (Not used by default)')
    parser.add_argument(
        '-Alpha', '--alpha_weight_regularization', 
        type=float, 
        required=False,
        default=0, 
        help='Experimental: Alpha for regularizing imputed C values. (default = 0, not used) ')

    return parser

def main():
    args = get_parser().parse_args(sys.argv[1:])

    C = np.load(args.reference, allow_pickle=True)
    _C = np.load(args.reference_partial, allow_pickle=True)
    B = np.load(args.bulk, allow_pickle=True)
    c = args.threads
    t = args.max_run_time_per_gurobi_instance

    # load weight_per_patient if it exists, else all ones
    if args.weight_per_patient:
        W = np.load(args.weight_per_patient, allow_pickle=True)
    else:
        W = np.ones(B.shape[0])

    # load alpha_weight_regularization if it exists, else 0
    if args.alpha_weight_regularization:
        Alpha = args.alpha_weight_regularization
    else:
        Alpha = 0

    C_initial = deepcopy(C)

    P_list = []
    C_list = []
    # for i in range(args.iters):

    i = 0
    fro_B_squared = np.square(LA.norm(B, 'fro'))
    delta = 1
    _lambda = np.inf
    while (delta > args.delta_percent_obj_error_threshold or _lambda > args.lambda_diff_successive_delta) \
            and (i < args.iters):

        P, _ = LP_get_P(C,B,c,t,W)
        P_list.append(P)
        C, obj_value = LP_get_C(P,B,_C, c,t,W,Alpha,C_initial)
        C_list.append(C)

        i += 1
        x = obj_value/fro_B_squared
        _lambda = delta - x
        delta = x
        print(f'main: After iteration = {i}')
        print(f'main: delta = {delta}')
        print(f'main: _lambda = {_lambda}')

    
    np.save(args.output_P,np.array(P_list, dtype=object), allow_pickle=True)
    np.save(args.output_C,np.array(C_list, dtype=object), allow_pickle=True)


if __name__=="__main__":
    #parser
    main()
