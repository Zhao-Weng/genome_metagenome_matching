import random


lines_of_file = 150366152
length = lines_of_file / 4

my_randoms = random.sample(range(0,length), 1000)
my_randoms = sorted(my_randoms)
print my_randoms


def parse(file_path):
	line_count = 0
	index = 0
	return_list = []

	with open(file_path) as fp:
	    for line in fp:

	    	if index == len(my_randoms):
	    		break

	    	line_count +=1
	    	if line_count == my_randoms[index] * 4 - 2:
	    		line = line.split()[0]
	    		return_list.append(line)
	    		index = index + 1

	print return_list


parse('TS_1_1.fastq')
