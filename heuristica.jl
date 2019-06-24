using Random
using Statistics

makespan_table = Dict{Array{Int64},Int64}()
last_line_table = Dict{Array{Int64},Array{Int64}}()

NUMBER_CANDIDATES = 270
STOP_GRASP = 1500
STOP_HILL_CLIMBING = 50

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

function makespan(sch::Array{Int64}, t, save_to_dict::Bool=true, use_last_line::Bool=false, previous_solution=nothing) # computes the makespan of a certain solution
    if haskey(makespan_table, sch) # if the dictionary already has the computed makespan
        return makespan_table[sch] # returns it and avoid computing it again
    
    # =========================================================================================== #    
    # este condicional não é mais utilizado pois adotamos um método alternativo apra construção da solução inicial
    # preservado para demonstração da implementação da otimização pensada
    elseif use_last_line == true && haskey(last_line_table, previous_solution) # when it intended to reuse the last line of the previous solution and that is already stored
        previous_last_line = last_line_table[previous_solution]                # retrieves the array corresponding to the last line of the previous solution
        current_last_line::Array{Int64} = [previous_last_line[1] + t[sch[length(sch)]],1] # current last line equivalent initialized as usual
                                                                                          # being the last job started on the first machine as soon as it finished on the previous (previous last line)
        number_machines = length(current_last_line)

        for m in 2:number_machines          # the last job at the current machine starts when both 
            top = previous_last_line[m-1]   # the same job has finished on the previous machine
            side = current_last_line[m-1]   # and the previous job has finished on the same machine
            current_last_line[m] = max(top,side) + current_last_line[m] 
        end

        last_line_table[sch] = current_last_line # stores the value for reuse

        return current_last_line[number_machines]
    # =========================================================================================== #

    else 
        number_jobs, number_machines = size(t)

        jobs_in_solution = length(sch)

        # t_line is the original t matrix with the order of the jobs alteres to represent the proposed schedule (sch)
        t_line = zeros(Int64, number_jobs, number_machines) 
        for i in 1:jobs_in_solution
            for j in 1:number_machines
                t_line[i,j] = t[sch[i],j] # puts the values in the positions referenced by the schedule
            end 
        end

        # finish_time is a matrix with the acumulated times of each job, showing exactly when job i in mahine j with finish, respecting the constrains
        finish_time = zeros(Int64, jobs_in_solution, number_machines)
        finish_time[1,1] = t_line[1,1] # only the job 1 at machine 1 never changes

        for i in 2:jobs_in_solution # fills each item of the first column with the sum of the values before in the t matrix
            finish_time[i,1] = t_line[i,1] + finish_time[i-1,1] # that is the computation of the finish times in the first machine, that will always execute in order
        end                                                     # starting a new job as soon as the previous finishes

        for j in 2:number_machines # fills the fist line with the sum of the previous values in the t matrix
            finish_time[1,j] = finish_time[1,j-1] + t_line[1,j] # that is the computation of the finish time of the first job of each machine
        end                                                     # always starting as soon as the first job of the previous machine finishes

        for i in 2:jobs_in_solution
            for j in 2:number_machines     
                side = finish_time[i,j-1] # each in between job depends on both the previous job in the same machine finishing and the same job in the previous machine finishg 
                top  = finish_time[i-1,j] # so its finish time will be the max value between those two
                finish_time[i,j] = max(top,side) + t_line[i,j] # + its own time to run
            end
        end

        makespan = finish_time[jobs_in_solution,number_machines] # the total makespan will be the finish time of the last job in the last machine

        if save_to_dict == true
            makespan_table[sch] = makespan # updates the makespan disctionary so that the recomputation can be avoided
        end

        if use_last_line == true
            last_line_table[sch] = finish_time[jobs_in_solution, :] # saves the last line of finish_time to be reused
        end

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

# ============================================================================================ #
#      Método orignal para a construção da colução inicial, não é mais utilizado               #
#                   código preservado para demonstração da implementação                       #

function compute_c_max(candidate::Array{Int64,1}, job_position::Int64, t::Array{Int64,2})
    number_jobs, number_machines = size(t)
    c_max::Int64 = 999999999

    if candidate == []           # if the candidate solution has no jobs yet
        return t[job_position,1] # returns its time in the first machine
    else
        temp_candidate = push!(copy(candidate),job_position) 

        return makespan(temp_candidate, t, false, true, candidate) # return the makespan of the jobs in the candidate + the job being tested as the schedule
    end

end

function randomized_greedy_construct(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    jobs_not_filled::Array{Int64,1} = [x for x in 1:number_jobs]
    candidate::Array{Int64,1} = []

    while length(jobs_not_filled) > 0
        execution_times = [(compute_c_max(candidate,i,t),i) for i in jobs_not_filled] # computes the makesan time for each job added to the current colution        
        execution_times = sort(execution_times)                     # orders those aditions by makespan
        next_job = ceil(Int,rand(1:(length(execution_times)*alpha+1))) # chooses a random job to be added among the best alpha%
        push!(candidate, execution_times[next_job][2])              # adds the job to the candidate being made
        jobs_not_filled = filter(x -> x != execution_times[next_job][2], jobs_not_filled) # removes job from the to be filled list 
    
    end

    return candidate
end
#                                                                                              #
# ============================================================================================ #

function modified_randomized_greedy_construct(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
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

    for k in 1:STOP_GRASP
        if k % 100 == 0
            println(k)
        end
        s_line = modified_randomized_greedy_construct(alpha, number_jobs, t) # gets the best solution given by the greedy construct

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

        println(s_star)
        println("initial solution makespan = ", makespan(initial_solution,t))
        println("s star makespan = ", makespan(s_star,t))
        println("execution time = ", execution_time)
    end
end

main()
