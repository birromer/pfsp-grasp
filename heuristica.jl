using Random

makespan_table = Dict{Array{Int64},Int64}()
NUMBER_CANDIDATES = 100
STOP_GRASP = 100
STOP_HILL_CLIMBING = 10


function read_instances(filename::String)
    sch = Array{Int64,1}

    stream = open(filename,"r")

    # splits first line by space and removes empy strings for better usage
    first_line = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
    number_jobs::Int64     = parse(Int64,first_line[1])
    number_machines::Int64 = parse(Int64,first_line[2]) 

    t   = zeros(Int64, number_jobs, number_machines)
    sch = [x for x in 1:number_jobs] # sets base schedule array with machine jobs in order

    # creation of the t matrix from the file
    for i in 1:number_jobs
        machine_values = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
        for j = 2:2:length(machine_values) # iterates jumping the index values
            k = Int(j/2)
            t[i,k] = parse(Int64,machine_values[j])
        end
    end

    return sch, t
end

function makespan(sch::Array{Int64}, t) # computes the makespan of a certain solution

    if haskey(makespan_table, sch) # if the dictionary already has the computed makespan
        return makespan_table[sch] # returns it and avoid computing it again
    else
        number_jobs, number_machines = size(t)

        # t_line is the original t matrix with the order of the jobs alteres to represent the proposed schedule (sch)
        t_line = zeros(Int64, number_jobs, number_machines) 
        for i in 1:number_jobs
            for j in 1:number_machines
                t_line[i,j] = t[sch[i],j] # puts the values in the positions referenced by the schedule
            end 
        end

        # finish_time is a matrix with the acumulated times of each job, showing exactly when job i in mahine j with finish, respecting the constrains
        finish_time = zeros(Int64, number_jobs, number_machines)
        finish_time[1,1] = t_line[1,1] # only the job 1 at machine 1 never changes

        for i in 2:number_jobs # fills each item of the first column with the sum of the values before in the t matrix
            finish_time[i,1] = t_line[i,1] + finish_time[i-1,1] # that is the computation of the finish times in the first machine, that will always execute in order
        end                                                     # starting a new job as soon as the previous finishes

        for j in 2:number_machines # fills the fist line with the sum of the previous values in the t matrix
            finish_time[1,j] = finish_time[1,j-1] + t_line[1,j] # that is the computation of the finish time of the first job of each machine
        end                                                     # always starting as soon as the first job of the previous machine finishes

        for i in 2:number_jobs
            for j in 2:number_machines     
                side = finish_time[i,j-1] # each in between job depends on both the previous job in the same machine finishing and the same job in the previous machine finishg 
                top  = finish_time[i-1,j] # so its finish time will be the max value between those two
                finish_time[i,j] = max(top,side) + t_line[i,j] # + its own time to run
            end
        end

        makespan = finish_time[number_jobs,number_machines] # the total makespan will be the finish time of the last job in the last machine

        makespan_table[sch] = makespan # updates the makespan disctionary so that the recomputation can be avoided

        return makespan
    end
end

function qsort!(array, lo, hi, t) # simple quicksort implementation to order candidate solutions in the randomized greedy construct
    i, j = copy(lo), copy(hi)
    while i < hi
        pivot = copy(array[(lo+hi)>>>1])
        while i <= j
            while makespan(array[i],t) < makespan(pivot,t); i = i+1; end
            while makespan(array[j],t) > makespan(pivot,t); j = j-1; end
            if i <= j
                array[i], array[j] = copy(array[j]), copy(array[i])
                i, j = i+1, j-1
            end
        end
        if lo < j; qsort!(array,lo,j,t); end
        lo, j = i, hi
    end
    return copy(array)
end

function generate_neighbour(solution::Array{Int64,1}) 
    number_jobs = length(solution) # the generation of a neighbour is the swap 
    i::Int32 = rand(1:number_jobs) # of two random jobs in the job order vector
    j::Int32 = rand(1:number_jobs)

    temp = solution[i]
    solution[i] = solution[j]
    solution[j] = temp

    return solution

end

function hill_climbing(solution::Array{Int64,1}, t) # hill climbing will be the local search used to improve the solution
    no_improvement_rounds::Int32 = 0 # counter for steps without improvement, will be used as stopping condition
    current_solution::Array{Int64,1} = deepcopy(solution) 
    current_makespan = makespan(solution,t) # precomputes the makespan for reuse 

    while no_improvement_rounds <= STOP_HILL_CLIMBING
        neighbour_solution = generate_neighbour(current_solution) # generates a new neighbour
        neighbour_makespan = makespan(neighbour_solution,t)       # and computes its makespan for reuse in case of a good one
    
        
        if neighbour_makespan < current_makespan # replaces the best current solution when has a shorter makespan
            current_makespan = neighbour_makespan
            current_solution = neighbour_solution
        else
            no_improvement_rounds += 1 # increment of the no improvement counter
        end
    end

    return current_solution
end

function initialize_candidate_set(number_jobs::Int64) # simple start for the candidate set
    candidates = Array{Int64}[]

    for i in 1:NUMBER_CANDIDATES
        new_candidate = randperm(number_jobs) # each candidate is a random permutation of the jobs
        push!(candidates, new_candidate)
    end

    return candidates
end


function randomized_greedy_construct(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    solutions = initialize_candidate_set(number_jobs) # initializes candidate set
    solutions = qsort!(solutions, 1, length(solutions), t) # orders set by makespan

    top_alpha = solutions[1:trunc(Int,alpha*length(solutions))] # takes the top alpha percentage of the results

    candidate = rand(top_alpha) # selcts random value between the top ones because sometimes the best one cant be much improved

    return candidate
end


function GRASP(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    s_star::Array{Int64,1} = randperm(number_jobs) # initialized optimal solution with random solution
    s_line = Array{Int64,1}

    for k in 1:STOP_GRASP
        
        #
        pre = makespan(s_star,t)
        #
        
        s_line = randomized_greedy_construct(alpha, number_jobs, t) # gets the best solution given by the greedy construct

        #
        pos = makespan(s_star,t)
        if (pre < pos)
            print("deu ruim")
        end
        println(makespan(s_star,t))
        #   

        s_line = hill_climbing(s_line,t) # tries to improve it with a local search



        if makespan(s_line,t) <= makespan(s_star,t) # if the new solution is better swaps it 
            println("troquei")
            s_star = copy(s_line)
        end

 

    end

    return s_star
end


function main()
    filename::String = ARGS[1]               # first parameter is the data file name 
    alpha::Float32 = parse(Float32, ARGS[2]) # second parameter is the alpha to be used in the construct
    Random.seed!(parse(Int64,ARGS[3]))       # third parameter is the randomness seed
    
    sch::Array{Int64,1}, t::Array{Int64,2} = read_instances(filename) # reads the instance from the file and loads initial schedule solution and t matrix

    number_jobs, number_machines = size(t)
    
    s_star = GRASP(alpha, number_jobs, t)

    println(s_star)
    println(makespan(s_star,t))
end

main()