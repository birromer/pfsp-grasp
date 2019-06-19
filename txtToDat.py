import os

filename = "VFR10_15_1_Gap"
file = open(filename + ".txt", "r")
nums = file.read()
i=0
j=0

# escrevendo N e M no arquivo:
while nums[i] != '\n':
    i += 1
    
nums_as_list = nums[j:i-1].split('  ')
j = i+1

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
for i in range(j, len(nums)):
    if nums[i] == '\n':
        nums_as_list = nums[j:i-1].split('  ')
        nums_as_list = [num for num in nums_as_list if num != '']
        nums_as_list_clean = [nums_as_list[i] for i in range(len(nums_as_list)) if i % 2 != 0]
        for k in range(1, m+1):
            if l == n and k == m:
                last_line = str(l) + " " + str(k) + " " + nums_as_list_clean[k-1] + ";\n\n\n"
            else:
                last_line = str(l) + " " + str(k) + " " + nums_as_list_clean[k-1] + "\n"
            f = open(filename + ".dat", "a")       
            f.write(last_line)
            f.close()
        l += 1
        j = i+1

f = open(filename + ".dat", "a")
f.write("end;")
f.close()