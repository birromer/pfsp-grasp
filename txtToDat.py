import os

filename = "VFR20_20_1_Gap"
file = open(filename + ".txt", "r")
nums = file.read()
i=0
j=0

# escrevendo N e M no arquivo:
while nums[i] != '\n':
    i += 1
    
nums_as_list = nums[j:i].split('  ')
g = i+1
j = g
f = open(filename + ".dat", "a")
f.write("data;\n\n")
f.write("param num_jobs := " + nums_as_list[0]+ ";\n\n")
f.write("param num_mach := " + nums_as_list[1] + ";\n\n")
f.write("param T := \n")
f.close()

n = int(nums_as_list[0])
m = int(nums_as_list[1])
k = 0
l = 1
for i in range(g, len(nums)):
    if nums[i] == '\n':
        nums_as_list = nums[j:i].split('  ')
        nums_as_list = [num for num in nums_as_list if num != '']
        nums_as_list_clean = [nums_as_list[i] for i in range(len(nums_as_list)) if i % 2 != 0]
        for k in range(m):
            if l == n and k == m-1:
                last_line = str(k+1) + " " + str(l) + " " + nums_as_list_clean[k] + ";\n\n\n"
            else:
                last_line = str(k+1) + " " + str(l) + " " + nums_as_list_clean[k] + "\n"
            f = open(filename + ".dat", "a")       
            f.write(last_line)
            f.close()
        l += 1
        j = i+1

f = open(filename + ".dat", "a")
f.write("end;")
f.close()
