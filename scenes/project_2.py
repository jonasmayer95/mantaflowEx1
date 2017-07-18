#
# Simple example scene for a 2D simulation
# Simulation of a buoyant smoke density plume with open boundaries at top & bottom
#
from manta import *

#this has to be done so keras doesn't kill itself
import locale
locale.setlocale(locale.LC_NUMERIC, "en_US.UTF-8")

import platform
import numpy as np
from keras.models import Sequential
from keras.layers import * 
from keras.models import load_model

#keras parameters
training = False

# solver params
compRes = 64
visRes = 128
compGs = vec3(compRes,compRes,1)
visGs = vec3(visRes,visRes,1)

sVis = Solver(name='avis', gridSize = visGs, dim=2)
s = Solver(name='comp', gridSize = compGs, dim=2)
s.timestep = 1
timings = Timings()

# prepare grids
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
density = s.create(RealGrid)
pressure = s.create(RealGrid)

visDensity = sVis.create(RealGrid)

bWidth=1
alpha=0.1
flags.initDomain(boundaryWidth=bWidth) 
flags.fillGrid()

setOpenBound(flags, bWidth,'yY',FlagOutflow|FlagEmpty) 

if (GUI):
	gui = Gui()
	gui.show( True ) 
	gui.pause()

source = s.create(Cylinder, center=compGs*vec3(0.5,0.1,0.5), radius=compRes*0.14, z=compGs*vec3(0, 0.02, 0))


def insert_dim(source, target, index):
	for x in range(source.shape[0]):
		for y in range(source.shape[1]):
			target[x][y][index] = source[x][y][0]
	return target;

def insert_dim_vec3(source, target, index):
	for x in range(source.shape[0]):
		for y in range(source.shape[1]):
			target[x][y][index] = source[x][y][0]
			target[x][y][index+1] = source[x][y][1]
			target[x][y][index+2] = source[x][y][2]
	return target;

	
#main loop
visArray= np.zeros((1,128,128,1))
compDensity = np.zeros((1,64,64,1))
compPressure = np.zeros((1,64,64,1))
compVel = np.zeros((1,64,64,3))

compArray= np.zeros((1,64,64,5))
model = load_model('smoke.h5')

for t in range(200):
	#mantaMsg('\nFrame %i' % (s.frame))

	if t<300:
		source.applyToGrid(grid=density, value=1)
		
	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2) 
	advectSemiLagrange(flags=flags, vel=vel, grid=vel,     order=2, openBounds=True, boundaryWidth=bWidth)
	
	resetOutflow(flags=flags,real=density) 

	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-2e-3,0), flags=flags)

	solvePressure(flags=flags, vel=vel, pressure=pressure)
	
	copyGridToArrayReal(density, compDensity[0])
	copyGridToArrayReal(pressure, compPressure[0])
	copyGridToArrayVec3(vel, compVel[0])

	compArray[0] = insert_dim(compDensity[0], compArray[0], 0)
	compArray[0] = insert_dim(compPressure[0], compArray[0], 1)
	compArray[0] = insert_dim_vec3(compVel[0], compArray[0], 2)
	
	visArray = model.predict(compArray, batch_size =1)
	copyArrayToGridReal(visArray[0], visDensity)
	
	#timings.display()    
	#gui.screenshot('plumeInstable\plume_%04d.png' % t) 
	s.step()
	sVis.step()

