
function read_instances(filename::String)
    sch = Array{Int64,1}

    stream = open(filename,"r")

    # splits first line by space and removes empy strings for better usage
    first_line = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
    number_jobs::Int64     = parse(Int64,first_line[1])
    number_machines::Int64 = parse(Int64,first_line[2]) 

    t   = zeros(Int64, number_jobs, number_machines)
    sch = [x for x in 1:number_machines] # sets base schedule array with machine numbers in order

    for i in 1:number_jobs
        machine_values = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
        for j = 2:2:length(machine_values) # iterates jumping the index values
            k = Int(j/2)
            t[i,k] = parse(Int64,machine_values[j])
        end
    end

    return sch, t
end


function makespan(sch::Array{Int64}, t_trans, makespan_table::Dict{Array{Int64},Int64}, alpha::Float32) # computes the makespan of a certain solution
    number_machines, number_jobs = size(t_trans)

    # t = t' # armazenar resultados da transposta

    t_line = zeros(Int64, number_machines, number_jobs)
    for i in 1:number_machines
        for j in 1:number_jobs
            t_line[i,j] = t_trans[sch[i],j]
        end 
    end

    finish_time = zeros(Int64, number_machines, number_jobs)
    finish_time[1,1] = t_line[1,1]

    for j in 2:number_jobs
        finish_time[1,j] = t_line[1,j] + finish_time[1,j-1]
    end

    for i in 2:number_machines
        finish_time[i,1] = finish_time[i-1,1] + t_line[i,1]
    end

    for i in 2:number_machines
        for j in 2:number_jobs
            top  = finish_time[i-1,j]
            side = finish_time[i,j-1]
            finish_time[i,j] = max(top,side) + t_line[i,j]
        end
    end

    makespan = finish_time[number_machines,number_jobs]

    makespan_table[sch] = makespan 

    return makespan
end


function hill_climbing(g::Function, solution, makespan_table::Dict{Array{Int64},Int64})

end


function initialize_candidate_set(number_solutions::Int32, number_jobs::Int64)
    candidates = Array{Array{Int64}}

    for i in 1:number_solutions:
        new_candidate = randperm(number_jobs)
        append(candidates, new_candidate)
    end

    return candidates
end


function randomized_greedy_construct(g::Function, makespan_table::Dict{Array{Int64},Int64}, alpha::Float32, number_candidates::int32, number_jobs::Int64)
    solutions = initialize_candidate_set(number_candidates, number_jobs)




end


function GRASP(alpha::Float32, stop::Int64)
    s_star::Intt64 = 0
    # s::Int64 = 0
    s_line::Int64 = 0

    for k in 1:stop
        s_line = construct(makespan, alpha, makespan_table)
        s_line = hill_climbing(makespan, s_line, makespan_table)
        if f(s_line) <= f(s_star)
            s_star = s_line
        end
    end
end

function main()
    filename::String = ARGS[1]
    alpha::Float32 = parse(Float32, ARGS[2])
    
    makespan_table = Dict{Array{Int64},Int64}()

    sch, t = read_instances(filename)
    
    t_trans = t'
    total_time = makespan(sch, t_trans, makespan_table, alpha)



    println("total time = ", total_time)

end

main()