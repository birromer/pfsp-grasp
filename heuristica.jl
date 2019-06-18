using DelimitedFiles



function read_instances(filename::String)
    sch = Array{Int32,1}
    
    stream = open(filename,"r")

    # splits first line by space and removes empy strings for better usage
    first_line = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
    number_machines::Int32 = parse(Int32,first_line[1])
    number_jobs::Int32 = parse(Int32,first_line[2]) 

    t = zeros(Int32, number_machines, number_jobs)

    sch = 1:number_machines  # sets base schedule array with machine numbers in order

    for i in 1:number_machines
        machine_values = filter((x) -> x != "", split(readline(stream, keep=false), ' '))
        for j = 2:2:length(machine_values) # iterates jumping the index items
            k = Int(j/2)
            t[i,k] = parse(Int32,machine_values[j])
        end
    end

    return sch, t
end

# function makespan(schedule, t)
#     time::Int32 = t[machine_order[1]][1]
#     for i in 2:length(schedule)
#         if 
#             time = time + 
#     end
# end


# function construct(g::Function, alpha::Float32)

# end


# function GRASP(alpha::Float32, stop::Int32)
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


end

main()