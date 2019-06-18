# use --tmlim nnn
param P := 32800;
param num_jobs;
param num_mach;
set M := 1..num_mach;
set N := 1..num_jobs;
param T{ r in M, i in N };

# tempo para completar a tarefa i na maquina r
var C{ r in M, i in N }  >= 0;
# 1 se tarefa i for feita antes da tarefa k, 0 cc
var D{ i in N, k in N } binary >= 0;
# makespan
var Cmax >= 0;

# funcao objetivo
minimize PFSP: Cmax;

s.t. executarTodas{ i in N }: C[1,i] >= T[1,i];
s.t. executarAteOFimEmOrdem{ r in M, i in N: r > 1 }: C[r,i] - C[r-1,i] >= T[r,i];
s.t. fixaOrdem1{ r in M, i in N, k in N: k > i and k > 1 }: C[r,i] - C[r,k] - P*D[i,k] >= T[r,i];
s.t. fixaOrdem2{ r in M, i in N, k in N: k > i and k > 1 }: C[r,i] - C[r,k] + P*D[i,k] <= P - T[r,k];
s.t. makespan{ r in M, i in N }: Cmax >= C[num_mach,i];

end;
