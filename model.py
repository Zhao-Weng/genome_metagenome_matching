import pandas
import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD
from keras.preprocessing import text
from keras.regularizers import l2, l1
from keras.models import model_from_json
from sklearn.metrics import precision_recall_fscore_support as score
import pdb

'''
input: truthTableName, preditedTableName, scoreFile
output: None
side effect: compare real label with predicted label and write preceision, recall, 
			 f1 score into file
'''
def bowtiePerformance(truthTableName, preditedTableName, scoreFile):
	truthFile = open(truthTableName, 'r')
	predictedFile = open(preditedTableName, 'r')
	truthArr = []
	predictedArr = []
	for truthLine, predictedLine in zip(truthFile, predictedFile):
		truthLine = truthLine[:-1] if truthLine[-1] == '\n' else truthLine
		predictedLine = predictedLine[:-1] if predictedLine[-1] == '\n' else predictedLine
		truthCurRow = truthLine.split('\t')
		predictedCurRow = predictedLine.split('\t')
		truthArr += (list(map(int, truthCurRow[1:])))
		predictedArr += (list(map(int, predictedCurRow[1:])))

	precision, recall, fscore, support = score(np.array(truthArr), np.array(predictedArr))
	scoreFile = open(scoreFile, 'a')
	scoreFile.write('{0}, {1}, {2}, {3}\n'.format(precision, recall, fscore, support))
	scoreFile.close()

'''
input: featureFile, modelName, truthTableFile, scoreFile, flag that indicates whether to go through training
output: None
side effect: read kmer features from featureFile， use ML model for training，save model, write to scoreFile
			 evaluation metrics
'''
def training(featureFile, modelName, truthTableFile, scoreFile, dotrain):
	dataframe = pandas.read_csv(featureFile, header=None)
	data = dataframe.values
	X = data[:,1:-1]
	Y = data[:, -1]
	if (dotrain):
		inputDim = X.astype(int).shape[1]  # num of columns
		# Building model
		model = Sequential()
		if modelName.endswith('nn'):
			model.add(Dense(inputDim//10, input_dim=inputDim, init='normal', activation='relu'))
			model.add(Dropout(0.05))
			model.add(Dense(inputDim//100, init='normal', activation='relu'))
			model.add(Dense(1, activation='sigmoid',W_regularizer=l2(0.01)))

		elif modelName.endswith('sigmoid'):
			model.add(Dense(1, activation='sigmoid',W_regularizer=l2(0.01), input_dim=inputDim))
		# sgd = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.0)
		# model.compile(loss='mae',
		#               optimizer=sgd,  
		#               metrics=["mae"]) 
		model.compile(optimizer='rmsprop', loss='binary_crossentropy')
		model.fit(X, Y, validation_split=0.33, epochs=50, batch_size=10)

		

		# serialize model to JSON
		model_json = model.to_json()
		with open(modelName + ".json", "w") as json_file:
		    json_file.write(model_json)
		# serialize weights to HDF5
		model.save_weights(modelName + ".h5")
		print("Saved model to disk")
		 

	json_file = open(modelName + ".json", 'r')
	loaded_model_json = json_file.read()
	json_file.close()
	model = model_from_json(loaded_model_json)
	model.load_weights(modelName + ".h5")
	print("Loaded model from disk")

	# Evaluating the prediction
	prediction = model.predict_classes(X).reshape([160])   # [160, 1] to [160]
	Y = Y.astype(int)									   # type object to int
	accuracy = sum(prediction == Y) / len(prediction)

	# write training accuracy to score file 
	# precision, recall, fscore, support = score(Y, prediction)
	# scoreFile = open(scoreFile, 'a')
	# scoreFile.write('{0}, {1}, {2}, {3}, {4}\n'.format(truthTableFile, precision, recall, fscore, support))
	# scoreFile.close()

	# write to file
	outputFile = open(truthTableFile, 'w')
	for genomeI in range(40):
		name = data[genomeI * 4, 0].split(':')[0]
		outputFile.write(name + ' ')
		string = ''
		for metaI in range(4):
			string += str(prediction[genomeI * 4 + metaI]) + '\t'
		outputFile.write(string[:-1] + '\n')
	outputFile.write('accuracy: {0}\n'.format(accuracy))
	outputFile.close()

# sort truthtable file with genome name as key
def sortTruthFile(truthTableName):
	ls = []
	with open(truthTableName, 'r') as truthFile:
		for line in truthFile:
			line = line[:-1] if line[-1] == '\n' else line
			curRow = line.split('\t')
			ls.append(curRow)
	ls = sorted(ls, key=lambda item:item[0])
	pdb.set_trace()
	truthFile.close()
	with open('challenge.txt', 'w') as truthFile:
		
		for curRow in ls:
			string = ''
			for item in curRow:
				string += (item + '\t')
			truthFile.write(string[:-1] + '\n')
	truthFile.close()



if __name__  == '__main__':
	# bowtiePerformance('strains2_training_truth.txt', 'predictedTruth_end.csv', 'bowtieScore.csv')

	kmerLens = [31, 17, 7]
	modelTypes = ['nn', 'sigmoid']
	for kmerLen in kmerLens:
		for modelType in modelTypes:
			training('features/features{0}train.csv'.format(kmerLen), 'models/model{0}_{1}'.format(kmerLen, modelType), 
				     'truthTables/truthTable{0}_{1}.txt'.format(kmerLen, modelType), 'modelScoreFile.csv', True)

	# modelTypes = ['sigmoid']
	# for modelType in modelTypes:
	# 	training('features/features17challenge.csv', 'models/model17_{0}'.format(modelType), 
	# 			     'truthTables/truthTable17_{0}_challenge.txt'.format(modelType), 'modelScoreFile.csv', False)
	#sortTruthFile('truthTables/truthTable17_sigmoid_challenge.txt')