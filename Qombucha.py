#!/usr/bin/env python



import gurobipy as gp, numpy as np, math
import sys, os, argparse


def LP_get_P(C,B,c,t):
    #n,m = S.shape
    
    #C is the cell type matrix C_i is for cel type i 
  
    n,m = C.shape # n cell type, m features
    r,_m = B.shape # r patients, m features
    
    
    print(r,n,m)
    
    # assert B and C have the same number of columns
    assert m ==_m
    
    model = gp.Model()
    model.Params.Threads = c
    model.Params.TimeLimit = t

    # initialize variables for inferred fractions for all patients
    P = np.empty((r,n),dtype=object)
    for i in range(r):
        for j in range(n):
            P[i,j]= model.addVar(vtype= gp.GRB.CONTINUOUS,name='var_p_{}_{}'.format(i,j))
 
    # set constraints for inferred fractions

    for i in range(r):
        for j in range(n):
            model.addConstr(P[i,j]>=0.0)
            model.addConstr(P[i,j]<=1.0)
        model.addConstr(gp.quicksum(P[i]) ==1.0,name = 'constr_p_{}'.format(i))
   

    # set objective

    #initialize B_hat and store all the entries in it
       # set objective

    B_hat = np.matmul(P,C)
    #residue = np.abs(B- B_hat).flatten()
    residue = B- B_hat
    residue = residue.flatten()
    R = residue
    
    


    Obj = gp.quicksum(r**2 for r in residue)
    model.setObjective(Obj, gp.GRB.MINIMIZE)

    model.optimize()

    P_arr = np.empty((r,n),dtype=float)
#    _C_arr = np.empty((n,m),dtype=float)
    
    #print(P)
    for i in range(r):
        for j in range(n):
            P_arr[i,j] = P[i,j].X

    return P_arr


def LP_get_C(P,B,_C, c,t):
    #n,m = S.shape
    
    #C is the cell type matrix C_i is for cel type i 
  
    r_,n = P.shape # n cell type, m features
    r,m = B.shape # r patients, m features
    
    
    print(r,n,m)
    
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
                C[i,j] = model.addVar(vtype= gp.GRB.CONTINUOUS,name='var_c_{}_{}'.format(i,j))
 
    for i in range(n):
        for j in range(m):
            if np.isnan(_C[i,j]):
                #model.addConstr(_C[i,j] in {0.0,0.5,1})
                #set constraints to inferred methylation has to be 0, 0.5 or 1
                
                model.addConstr(C[i,j]>=0,name ='constr_c_{i}_{j}_lb')
                model.addConstr(C[i,j]<=1,name ='constr_c_{i}_{j}_ub')

    


    # set objective

    #initialize B_hat and store all the entries in it
       # set objective

    B_hat = np.matmul(P,C)
    #residue = np.abs(B- B_hat).flatten()
    residue = B- B_hat
    residue = residue.flatten()
    R = residue
    
    


    Obj = gp.quicksum(r**2 for r in residue)
    model.setObjective(Obj, gp.GRB.MINIMIZE)

    model.optimize()

    C_arr = np.empty((n,m),dtype=float)
    
    for i in range(n):
        for j in range(m):
            if np.isnan(_C[i,j]):
                C_arr[i,j] = C[i,j].X
            else:
                C_arr[i,j] = _C[i,j]
    
    return C_arr





def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-C', '--reference', type=str, required=True)
    parser.add_argument('-Cp', '--reference_partial', type=str, required=True)
    parser.add_argument('-B', '--bulk', type=str, required=True)
    parser.add_argument('-op', '--output_P', type=str, required=True)
    parser.add_argument('-oc', '--output_C', type=str, required=True)
 #   parser.add_argument('-q', '--q', type=float, required=True)
    parser.add_argument('-I', '--iters', type=int, required=False,default = 10)
    parser.add_argument('-c', '--threads', type=int, required=False,default=2)
    parser.add_argument('-t', '--run_time', type=int, required=False,default=100000)

    return parser

def main():
    args = get_parser().parse_args(sys.argv[1:])
 #   assert args.q <= 1

    C = np.load(args.reference, allow_pickle=True)['m']
    _C = np.load(args.reference_partial, allow_pickle=True)['m']
    B = np.load(args.bulk, allow_pickle=True)['m']
    c = args.threads
    t = args.run_time

    P_list = []
    C_list = []
    for i in range(args.iters):

        P = LP_get_P(C,B,c,t)
        P_list.append(P)
        C = LP_get_C(P,B,_C, c,t)
        C_list.append(C)
    
    #print(P)



    
    np.save(args.output_P,np.array(P_list, dtype=object), allow_pickle=True)
    np.save(args.output_C,np.array(C_list, dtype=object), allow_pickle=True)


if __name__=="__main__":
    #parser
    main()
