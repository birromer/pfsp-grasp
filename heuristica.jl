using Random
using Statistics

makespan_table = Dict{Array{Int64},Int64}()
NUMBER_CANDIDATES = 1000
STOP_GRASP = 1000
STOP_HILL_CLIMBING = 100

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
    i, j = lo, hi
    while i < hi
        pivot = array[(lo+hi)>>>1]
        while i <= j
            while makespan(array[i],t) < makespan(pivot,t); i = i+1; end
            while makespan(array[j],t) > makespan(pivot,t); j = j-1; end
            if i <= j
                array[i], array[j] = array[j], array[i]
                i, j = i+1, j-1
            end
        end
        if lo < j; qsort!(array,lo,j,t); end
        lo, j = i, hi
    end
    return array
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
    current_solution::Array{Int64,1} = copy(solution) 
    current_makespan = makespan(solution,t) # precomputes the makespan for reuse 

    while no_improvement_rounds <= STOP_HILL_CLIMBING
        neighbour_solution = generate_neighbour(copy(current_solution)) # generates a new neighbour
        neighbour_makespan = makespan(copy(neighbour_solution),t)       # and computes its makespan for reuse in case of a good one
        
        if neighbour_makespan < current_makespan # replaces the best current solution when has a shorter makespan
            current_makespan = copy(neighbour_makespan)
            current_solution = copy(neighbour_solution)
        else
            no_improvement_rounds += 1 # increment of the no improvement counter
        end
    end

    return current_solution
end

function initialize_candidate_set(number_jobs::Int64) # simple start for the candidate set
    p_solutions = Array{Int64}[]

    for i in 1:NUMBER_CANDIDATES
        new_candidate = randperm(number_jobs) # each candidate is a random permutation of the jobs
        push!(p_solutions, new_candidate)
    end

    return p_solutions
end


function randomized_greedy_construct(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    solutions = initialize_candidate_set(number_jobs) # initializes candidate set
    solutions = qsort!(copy(solutions), 1, length(solutions), t) # orders set by makespan

    top_alpha = solutions[1:trunc(Int,alpha*length(solutions))] # takes the top alpha percentage of the results

    candidate = rand(top_alpha) # selcts random value between the top ones because sometimes the best one cant be much improved

    return candidate
end


function GRASP(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    initial_solution::Array{Int64,1} = randperm(number_jobs) # random solution for initial one
    s_star::Array{Int64,1} = copy(initial_solution)          # initialized optimal solution with random solution
    s_line = Array{Int64,1}
    pre = 0
    pre_1 = randperm(number_jobs)

    for k in 1:STOP_GRASP
        s_line = randomized_greedy_construct(alpha, number_jobs, t) # gets the best solution given by the greedy construct

        s_line = hill_climbing(copy(s_line),t) # tries to improve it with a local search

        if makespan(s_line,t) < makespan(s_star,t) # if the new solution is better swaps it 
            s_star = copy(s_line)
        end
    end

    return initial_solution, s_star
end


function main()
    filename::String = ARGS[1]               # first parameter is the data file name 
    alpha::Float32 = parse(Float32, ARGS[2]) # second parameter is the alpha to be used in the construct
    Random.seed!(parse(Int64,ARGS[3]))       # third parameter is the randomness seed
    
    s_stars           = Array{Int64}[]
    initial_solutions = Array{Int64}[]

    time_elapsed              ::Array{Float64} = [] 
    initial_solutions_makespan::Array{Float64} = []
    s_stars_makespan          ::Array{Float64} = []
    
    sch::Array{Int64,1}, t::Array{Int64,2} = read_instances(filename) # reads the instance from the file and loads initial schedule solution and t matrix

    number_jobs, number_machines = size(t)
    
    for i in 1:1
        Random.seed!(parse(Int64,ARGS[3])*i)
        execution_time  = @elapsed initial_solution, s_star  = GRASP(alpha, number_jobs, t)
        
        push!(time_elapsed              , execution_time)
        push!(s_stars                   , s_star)
        push!(s_stars_makespan          , makespan(s_star,t))
        push!(initial_solutions         , initial_solution)
        push!(initial_solutions_makespan, makespan(initial_solution,t))
    end

    mean_initial_solution = mean(initial_solutions_makespan)
    mean_s_star           = mean(s_stars_makespan)
    std_dev_s_star        = std(s_stars_makespan)
    mean_execution_time   = mean(time_elapsed)

    println("mean initial solution = ", mean_initial_solution)
    println("mean best solution = ", mean_s_star)
    println("std dev best solution = ", std_dev_s_star)
    println("mean execution time = ", mean_execution_time, "\n")
end

main()