def parse(file_path):
	string = ''
	with open(file_path) as fp:
	    for line in fp:
	    	if line[0] == '>':
	    		continue
	    	else:
	    		line = line.split()[0]
	    		string = string + str(line)
	    	
	print string
     

parse('Alistipes_T_3.fna')

