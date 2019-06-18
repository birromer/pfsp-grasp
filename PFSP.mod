# sei la ta errado
param P := 32800;
param num_jobs;
param num_mach;
set N := 1..num_jobs;
set M := 1..num_mach;
param T{ r in M, i in N };

# tempo para completar a tarefa i na maquina r
var C{ r in M, i in N }  >= 0;
# 1 se tarefa i for feita antes da tarefa k, 0 cc
var D{ i in N, k in N } binary >= 0;
# makespan
var Cmax integer >= 0;

# funcao objetivo
maximize PFSP: -Cmax;

s.t. executarTodas{ r in M, i in N: i == 1 }: -C[r,i] <= -T[r,i];
s.t. executarAteOFimEmOrdem{ r in M: r > 1, i in N }: C[r-1,i] - C[r,i] <= -T[r,i];
s.t. naosei1{ r in M, i in N: i < num_jobs, k in N: k > i }: C[r,k] - C[r,i] - P*D[i,k] <= -T[r,i]
s.t. naosei2{ r in M, i in N: i < num_jobs, k in N: k > i }: C[r,i] - C[r,k] + P*D[i,k] <= P - T[r,k]
s.t. makespan{ r in M, i in N }: C[r,i] <= Cmax;

end;
