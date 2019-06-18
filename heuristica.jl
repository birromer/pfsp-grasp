
function read_instances(filename::String)
    sch = Array{UInt64,1}

    stream = open(filename,"r")

    # splits first line by space and removes empy strings for better usage
    first_line = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
    number_machines::UInt64 = parse(UInt64,first_line[1])
    number_jobs::UInt64 = parse(UInt64,first_line[2]) 

    t = zeros(UInt64, number_machines, number_jobs)

    sch = [x for x in 1:number_machines]#1:number_machines  # sets base schedule array with machine numbers in order

    for i in 1:number_machines
        machine_values = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
        for j = 2:2:length(machine_values) # iterates jumping the index items
            k = Int(j/2)
            t[i,k] = parse(UInt64,machine_values[j])
        end
    end

    return sch, t
end


function makespan(sch, t) # computes the makespan of a certain solution
    time_used::UInt64 = t[sch[1],1]
    number_machines, number_jobs = size(t)

    t_line = zeros(UInt64, number_machines, number_jobs)
    for i in 1:number_machines
        for j in 1:number_jobs
            t_line[i,j] = t[sch[i],j]
        end 
    end

    finish_time = zeros(UInt64, number_machines, number_jobs)

    for j in 1:number_jobs
        finish_time[1,j] = t_line[1,j]
    end

    for i in 2:number_machines
        finish_time[i,1] = finish_time[i-1,1] + t_line[i,1]
    end

    for i in 2:number_machines
        for j in 2:number_jobs
            top  = finish_time[i-1,j-1] + finish_time[i,j-1]
            side = finish_time[i-1,j-1] + finish_time[i-1,j]
            finish_time[i,j] = max(top,side)
        end
    end

    # for i = 1:number_machines
    #     for j = 1:number_jobs
    #         print(finish_time[i,j])
    #         print(' ')
    #     end
    #     print('\n')
    # end

    return finish_time[number_machines,number_jobs]
end



# function construct(g::Function, alpha::Float32)

# end


# function GRASP(alpha::Float32, stop::UInt64)
#     s_star::Float64 = 0
#     s::Float64 = 0
#     s_line::Float64 = 0

#     for k in 1:stop
#         s = construct(f, alpha)
#         s_line = local_search(f,s)
#         if f(s_line) <= f(s_star)
#             s_star = s_line
#         end
#     end
# end

function main()
    filename::String = ARGS[1]
    sch, t = read_instances(filename)
    total_time = makespan(sch, t)

    println("total time = ", total_time)

end

main()