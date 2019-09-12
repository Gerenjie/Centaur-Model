import io
import sys
import numpy as np

def getChebsFromFiles ():
	obs = []
	for fileNum in range(10):
		with open('chebs'+str(fileNum)+'.txt', 'r') as file:
			str=file.read()
			str=str.replace('\n', '')
			str=str.replace(' ' , '')
			str=str.replace('\t', '')
			abc = str.split('][')
			abc[0] = abc[0][1:]
			abc[-1] = abc[-1][:-1]
			for particle in abc:
				piecewiseChebStrings = particle.split('}{')
				piecewiseChebStrings[0] = piecewiseChebStrings[0] [1:]
				piecewiseChebStrings[-1]= piecewiseChebStrings[-1][:-1]
				obj = []
				for coordinatePiecewiseCheb in piecewiseChebStrings:
					pieces = coordinatePiecewiseCheb.split('),') 
					pieces[-1]= pieces[-1][:-1]
					newDict = {}
					for piece in pieces:
						key, value = piece.split(':')
						key = float(key)
						value = value[7:-1]
						coefficients = value.split(',')
						for i in range(len(coefficients)):
							coefficients[i] = float(coefficients[i])
						coefficients = np.array(coefficients)
						newDict[key]=coefficients
					obj.append(newDict)
				obs.append(obj)		
			print(float(obs[0][0][0.0][0]),'\n')
			print(type(float(obs[0][0][0.0][0])))
			print(len(obs))
		file.close()	
	return obs
