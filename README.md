# Scratch code for LMFDB

download_database.py queries a database of generating vectors for all generating vectors in a certain range by genus and downloads them into 10 separate files, split as evenly as possible, with the necessary information for the no-dp.py file. 

no-db.py performs the magma calculations for braid and topological equivalence without interacting with the Mongo databse. It takes a single command-line argument specifying an input file which is expected to be the output from download_database.py. The code will then output its results to an output file named the same as the input file but ending in -output.yml. 

upload-orbit.py uploads the output files of no-dp.py to the database of generating vectors. It takes the names of files as command-line arguments, and can take several at once.
