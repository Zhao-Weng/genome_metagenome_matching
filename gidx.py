import random, os, csv, pdb

def parseMeta(file_path, num):
  """
  a function to parse the metagenome. Extract the sequence string and return a list of strings.

  param:
    file_path: path of the metagenome file
    num: the number of metagenome file
  """
  if num == 1:
    lines_of_file = 150366152
  elif num == 2:
    lines_of_file = 174516324
  elif num == 3:
    lines_of_file = 266426456
  elif num == 4:
    lines_of_file = 432415648 
  length = lines_of_file / 4

  my_randoms = random.sample(range(0,length), 500)
  my_randoms = sorted(my_randoms)

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

  return return_list

def parseGe(file_path):
  """
  a function to parse the genome. Get rid of the node id and concatenate all the sequences into 
  a long string. Return the concatenated sequence string

  param:
    file_path: path of the genome file
  """
  string = ''
  with open(file_path) as fp:
    for line in fp:
      if line[0] == '>':
        continue
      else:
        line = line.split()[0]
        string = string + str(line)
  return string

def buildIndex(genome, kmer):
  index = set()
  for i in range(0, len(genome)-kmer+1):
    index.add(genome[i:i+kmer])
  return index

def buildDatapoint(genome, mglist, kmer=31):
  """
  a function to build a feature vector for a single genome vs a single metagenome

  param:
    genome: a concatenated string
    mglist: a list of metagenome reads
    kmer=31, default value
  """
  ## build hash index##
  index = buildIndex(genome, kmer)
  
  lengths = [len(s) for s in mglist]
  n_kmer = max(lengths)-kmer+1
  n_reads = len(mglist)
  matrix = [[0 for i in range(n_kmer)] for j in range(n_reads)]
  for r in range(n_reads):
    read = mglist[r]
    for i in range(n_kmer):
      if i >= len(read):
        break
      if read[i:i+kmer] in index:
        matrix[r][i] = 1
  dp = []
  for line in matrix:
    dp += line
  return dp

def readTruthtable(truth_table):
  with open(truth_table, 'r') as f:
    tt = f.readlines()
  genomes = ['' for i in range(40)]
  meta_genomes = ['TS_1', 'TS_2', 'TS_3', 'TS_4']
  truth_mat = [[0, 0, 0, 0] for i in range(40)]
  for i in range(len(tt)):
    line = tt[i].strip('\n').split('\t')
    genomes[i] = line[0]
    truth_mat[i][0] = int(line[1])
    truth_mat[i][1] = int(line[2])
    truth_mat[i][2] = int(line[3])
    truth_mat[i][3] = int(line[4])
  return genomes, meta_genomes, truth_mat

def writeCSV(genomes_path, metagenomes_path, truth_path):
  """
  a function to generate all the feature vectors for each genome vs each metagenome.
  It will write all the features to a csv file. The csv file will have 160 rows, with 
  each row as a data point. For each row, it will be made up as three part:
    First entry: genome_name:metagenome_name
    Last entry: groundtruth label
    Rest: feature values

  genomes_path: the path of the folder that contains all the genomes
  metagenomes_path: the path of the folder that contains all the metagenomes
  truth_path: the path of the groundtruth table
  """
  genomes_names, mg_names, truth_mat = readTruthtable(truth_path)
  with open('features_17mer.csv', 'w') as f:
    wr = csv.writer(f)

    for i in range(len(genomes_names)):
      genome_name = genomes_names[i]
      genome_path = os.path.join(genomes_path, genome_name+'.fna')
      print('Parsing '+genome_name+'...')
      genome = parseGe(genome_path)

      for j in range(len(mg_names)):
        mglist = []
        if os.path.isfile(mg_names[j]):
          with open(mg_names[j], 'rb') as pickle_file:
            mglist = pickle.load(pickle_file)
        else:
          mg_name1 = mg_names[j]+'_1'
          mg_name2 = mg_names[j]+'_2'
          mg_path1 = os.path.join(metagenomes_path, mg_name1+'.fastq')
          mg_path2 = os.path.join(metagenomes_path, mg_name2+'.fastq')

          print('Parsing '+mg_names[j]+'...')
          mglist1 = parseMeta(mg_path1, j+1)
          mglist2 = parseMeta(mg_path2, j+1)
          mglist = mglist1+mglist2
          with open(mg_names[j], 'wb') as pickle_file:
            pickle.dump(mglist, pickle_file)

        print('Building feature '+genome_name+':'+mg_names[j]+'...')
        feature = buildDatapoint(genome, mglist, 17)
        dp = [genome_name+':'+mg_names[j]]+feature+[truth_mat[i][j]]
        wr.writerow(dp)

writeCSV('strains2_training_genomes', 'strains2_training_metagenomes', 'strains2_TRAINING_truth.txt')