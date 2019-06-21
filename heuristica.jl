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
    sch = [x for x in 1:number_jobs] # sets base schedule array with machine numbers in order

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

    if haskey(makespan_table, sch)
        return makespan_table[sch]
    else
        number_jobs, number_machines = size(t)

        t_line = zeros(Int64, number_jobs, number_machines)
        for i in 1:number_jobs
            for j in 1:number_machines
                t_line[i,j] = t[sch[i],j]
            end 
        end

        finish_time = zeros(Int64, number_jobs, number_machines)
        finish_time[1,1] = t_line[1,1]

        for i in 2:number_jobs
            finish_time[i,1] = t_line[i,1] + finish_time[i-1,1]
        end

        for j in 2:number_machines
            finish_time[1,j] = finish_time[1,j-1] + t_line[1,j]
        end

        for i in 2:number_jobs
            for j in 2:number_machines
                side = finish_time[i,j-1]
                top  = finish_time[i-1,j]
                finish_time[i,j] = max(top,side) + t_line[i,j]
            end
        end

        makespan = finish_time[number_jobs,number_machines]

        makespan_table[sch] = makespan 

        return makespan
    end
end

function qsort!(array, lo, hi, t)
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
    return array
end

function generate_neighbour(solution::Array{Int64,1})
    number_jobs = length(solution)
    i::Int32 = rand(1:number_jobs)
    j::Int32 = rand(1:number_jobs)

    temp = solution[i]
    solution[i] = solution[j]
    solution[j] = temp

    return solution

end

function hill_climbing(solution::Array{Int64,1}, t)
    no_improvement_rounds::Int32 = 0
    current_solution::Array{Int64,1} = copy(solution)
    current_makespan = makespan(solution,t)

    while no_improvement_rounds <= STOP_HILL_CLIMBING
        neighbour_solution = generate_neighbour(current_solution)
        neighbour_makespan = makespan(neighbour_solution,t)
    
        
        if neighbour_makespan < current_makespan
            current_makespan = neighbour_makespan
            current_solution = neighbour_solution
        else
            no_improvement_rounds += 1
        end
    end

    # println(current_solution)
    # println("current_makespan = ", current_makespan)
    # println("makespan(current_solution = ", makespan(current_solution,t))

    return current_solution
end

function initialize_candidate_set(number_jobs::Int64)
    candidates = Array{Int64}[]

    for i in 1:NUMBER_CANDIDATES
        new_candidate = randperm(number_jobs)
        push!(candidates, new_candidate)
    end

    return candidates
end


function randomized_greedy_construct(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    solutions = initialize_candidate_set(number_jobs)
    solutions = qsort!(solutions, 1, length(solutions), t)

    top_alpha = solutions[1:trunc(Int,alpha*length(solutions))]

    candidate = rand(top_alpha)

    return candidate
end


function GRASP(alpha::Float32, number_jobs::Int64, t::Array{Int64,2})
    s_star::Array{Int64,1} = randperm(number_jobs)
    s_line = Array{Int64,1}

    
    for k in 1:STOP_GRASP
        pre = makespan(s_star,t)
        s_line = randomized_greedy_construct(alpha, number_jobs, t)

        s_line = hill_climbing(s_line,t)



        if makespan(s_line,t) <= makespan(s_star,t)
            println("troquei")
            s_star = copy(s_line)
        end
        pos = makespan(s_star,t)
        if (pre < pos)
            print("deu ruim")
        end

        println(makespan(s_star,t))
    end

    return s_star
end


function main()
    Random.seed!(parse(Int64,ARGS[3]))
    filename::String = ARGS[1]
    alpha::Float32 = parse(Float32, ARGS[2])
    
    sch::Array{Int64,1}, t::Array{Int64,2} = read_instances(filename)

    number_jobs, number_machines = size(t)
    
    s_star = GRASP(alpha, number_jobs, t)

    println(s_star)
    println(makespan(s_star,t))
end

main()