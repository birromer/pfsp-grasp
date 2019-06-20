using Random

makespan_table = Dict{Array{Int64},Int64}()
NUMBER_CANDIDATES = 100
STOP = 100

function read_instances(filename::String)
    sch = Array{Int64,1}

    stream = open(filename,"r")

    # splits first line by space and removes empy strings for better usage
    first_line = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
    number_jobs::Int64     = parse(Int64,first_line[1])
    number_machines::Int64 = parse(Int64,first_line[2]) 

    t   = zeros(Int64, number_jobs, number_machines)
    sch = [x for x in 1:number_jobs] # sets base schedule array with machine numbers in order
    # temp = sch[1]
    # sch[1] = sch[3]
    # sch[3] = temp 

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

    for i in 1:number_jobs
        for j in 1:number_machines
            print(finish_time[i,j], ' ')
        end
        print('\n')
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

function qsort!(a,lo,hi)
    i, j = lo, hi
    while i < hi
        pivot = a[(lo+hi)>>>1]
        while i <= j
            while makespan(a[i]) < makespan(pivot); i = i+1; end
            while makespan(a[j]) > makespan(pivot); j = j-1; end
            if i <= j
                a[i], a[j] = a[j], a[i]
                i, j = i+1, j-1
            end
        end
        if lo < j; qsort!(a,lo,j); end
        lo, j = i, hi
    end
    return a
end

function hill_climbing(g::Function, solution)

end

function initialize_candidate_set(number_jobs::Int64)
    candidates = Array{Int64}[]

    for i in 1:NUMBER_CANDIDATES
        new_candidate = randperm(number_jobs)
        push!(candidates, new_candidate)
    end

    for i in candidates
        println(i)
    end

    return candidates
end


function randomized_greedy_construct(g::Function, alpha::Float32, number_jobs::Int64)
    solutions = initialize_candidate_set(number_jobs)
    solutions = qsort!(solutions, 1, length(solutions))

    return candidate
end


function GRASP(alpha::Float32, number_jobs)
    s_star::Int64 = 0
    s_line::Int64 = 0

    for k in 1:STOP
        s_line = randomized_greedy_construct(makespan, alpha, number_jobs)
        s_line = hill_climbing(makespan, s_line, makespan_table)
        if f(s_line) <= f(s_star)
            s_star = s_line
        end
    end
end

function main()
    filename::String = ARGS[1]
    alpha::Float32 = parse(Float32, ARGS[2])
    
    sch::Array{Int64,1}, t::Array{Int64,2} = read_instances(filename)

    number_jobs, number_machines = size(t)
    
    # total_time = makespan(sch, t)

    s_star = GRASP(alpha, number_jobs)


end

main()