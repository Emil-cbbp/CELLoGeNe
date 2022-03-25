import lgnReader
import Network
import sys
import time


""" Prints the time elapsed since the given start time. """
def time_elapsed(start_time):
	run_time = time.time() - start_time
	hours = run_time // 3600
	minutes = (run_time % 3600) // 60
	seconds = run_time - hours * 3600 - minutes * 60
	print("--- Time elapsed: %d h %d min and %.3f s --- Current time: %s ---"
		% (int(hours), int(minutes), seconds, time.strftime("%H:%M:%S", time.localtime())))



"""
	The main method of the program. Creates a reader and reads a lgn-file,
	sets up the network, and starts the run
"""
if(__name__ == "__main__"):
	start = time.time()
	
	""" If no argument (.lgn file) provided while running, use the demo file. """
	if(len(sys.argv) == 1):
		read = lgnReader.Reader("HowTo.lgn")
	else:
		read = lgnReader.Reader(sys.argv[1])
	geneData = read.geneData
	commands = read.commands
	lgn = Network.Network(geneData, commands)
	if(commands['runType'] == 'E'):
		lgn.runExhaustive()
	elif(commands['runType'] == 'M'):
		if('frac' in commands):
			lgn.runMonteCarlo(frac=commands['frac'])
		elif('tot' in commands):
			lgn.runMonteCarlo(tot=commands['tot'])
		else:
			print('Keyword is missing')

	time_elapsed(start)
