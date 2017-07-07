Here's a quick guide how to run this example. Note that the scripts here assume
that the executable is located in ".../mantaflow/build/", and most scripts were
primarily developed for Unix systems (Windows will work, but require path modifications)..

- first generate data by calling: 
    >>> manta manta_genSimData.py
  (windows, e.g.: ../../build/Release/manta.exe manta_genSimData.py)

- you can also use the tf_genManySims.py script to generate 10 data sets in one go. Note, manta_X.py files
  are intended to be run with mantaflow, while tf_Y.py files should be run with python, i.e. use
  "python tf_genManySims.py" to generate the data.

- then use it to train a first model with tensorflow. This will take ca. 2 min. E.g., with: 
    >>>python tf_train_pressure.py out 0 fromSim 1001 toSim -1 useVelocities 1 outName pressure bWidth 1 trainingEpochs 50000
  now you should have a trained model checkpoint for test_0001 in the data directory (../data by default). 
  outName denotes the name of the learned quantity. The input data name should match.

- you can then use the model to generate an output with the command below. It assumes the trained
  model was generated in ../data/test_0001 , and that the sim data to apply it to is number 1000:
    >>>python tf_train_pressure.py out 1 fromSim 1000 toSim -1 useVelocities 1 outName pressure bWidth 1 loadModelTest 1 brightenOutput 10
