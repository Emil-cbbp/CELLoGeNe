import numpy as np
import matplotlib.pyplot as plt

OE, KD, FREE, M, A, ON, OFF = +1, -1, 0, +2, -2, +3, -3

class Cell:
	def __init__(self, init_state, genes, gene_statuses, landscape, dim, beta=1):
		self.s_init = init_state#int(init_state, 2)
		self.genes = np.array(genes)
		self.gene_statuses = gene_statuses
		self.L = np.array(landscape)
		self.dim = dim
		self.s = self.s_init
		self.s_history = [self.s_init]
		self.beta = beta
	
	
	
	
	def energy(self, s):
		if(isinstance(s, (int, np.integer))):
			return self.L[s]
		elif(isinstance(s, str)):
			s = int(s, 2)
			return self.L[s]
		elif(isinstance(s, list)):
			states = s
			for i, s in enumerate(states):
				if(isinstance(s, (int, np.integer))):
					continue
				elif(isinstance(s, str)):
					states[i] = int(s, 2)
				else:
					print(("The cell state is not given in a valid format. " +
					"Must be int, str, not {:s}").format(str(type(s))))
			states = np.array(states)
			return self.L[states]
		
		elif(isinstance(s, np.ndarray)):
			states = s
			if(isinstance(states[0], (int, np.integer))):
				return(self.L[states])
			elif(isinstance(states[0], str)):
				states = np.array([int(x, 2) for x in states])
				return(self.L[states])
			else:
				print('Here?')
				print(("The cell state is not given in a valid format. " +
				"Must be int, str, not {:s}").format(str(type(states[0]))))
		else:
			print(("The cell state is not given in a valid format. " +
					"Must be int, str, or list/np.array of those, not {:s}").format(str(type(s))))
	
	
	def find_neighbours(self, s, asnp=False):
		if(isinstance(s, (int, np.integer))):
			pass
		elif(isinstance(s, str)):
			s = int(s, 2)
		neigh = [s ^ (1 << m) for m in range(self.dim)]
		if(asnp):
			return np.array(neigh, dtype=int)
		else:
			return neigh
	
	def take_step(self, no_stay=False):
		neigh = self.find_neighbours(self.s)
		if(no_stay):
			inc_gs = -1
		else:
			neigh.append(self.s)	# Add the state to list of neighbours to enable to stay at same state
			inc_gs = len(self.gene_statuses)
		
		E_2 = self.energy(neigh)
		E_1 = self.energy(self.s)
		
		
		if(no_stay):
			a = np.exp(-(E_2 - E_1) * self.beta) * self.gene_statuses[:-1]
		else:
			a = np.exp(-(E_2 - E_1) * self.beta) * self.gene_statuses
		
		
		a_0 = np.sum(a)
		r = np.random.rand()
		a_cum = np.cumsum(a)
		
		mu = np.min(np.where(a_cum >= (r * a_0)))
		
		self.s = neigh[mu]
		self.s_history.append(self.s)
	
	def simulate_cell(self, stop=-1, put_lim=3, perturb_init=False):
		put = 0
		prev_s = self.s_init
		steps = 0
		if(perturb_init):
			self.perturb()
			prev_s = self.s
		while(put < put_lim):
			self.take_step()
			if(self.s == prev_s):
				put += 1
			else:
				put = 0
			prev_s = self.s
			steps += 1
			if(stop >= 0 and steps >= stop):
				break
	
	def perturb(self, no_stay=True):
		if(no_stay):
			self.take_step(no_stay=no_stay)
		
		
		
	

class Simulation:
	def __init__(self, N, s_init, genes, gene_statuses, landscape, beta=1):
		self.N = N
		self.s_init = int(s_init, 2)
		self.genes = np.array(genes)
		gene_statuses.append(1)
		self.gene_statuses = np.array(gene_statuses)
		self.L = np.array(landscape)
		self.dim = len(genes)
		self.beta = beta
	
		tmp = bin(len(self.L))[2:]
		tmp = len(tmp) - len(tmp.rstrip('0'))
		if(tmp != self.dim):
			print("Number of genes and dimension of landscape are not matching.")

		self.cells = [Cell(self.s_init, self.genes, self.gene_statuses, self.L, self.dim, self.beta) for i in range(self.N)]
	
	def simulate_cells(self, stop=-1, perturb_init=False):
		for c in self.cells:
			c.simulate_cell(stop, perturb_init=perturb_init)
	
	def end_states(self):
		s_fin = []
		for i, c in enumerate(self.cells):
			s_fin.append(c.s_history[-1])
		uni = np.unique(s_fin, return_counts=True)
		tmp = np.zeros((len(uni[0]), 2), dtype=int)
		tmp[:,0] = uni[0]
		tmp[:,1] = uni[1]
		self.uni = tmp
		self.uni = self.uni[np.argsort(self.uni[:,1])][::-1]
		self.uni_dic = {np.binary_repr(d[0], self.dim): d[1] / self.N for d in self.uni}
	
	def get_end_states(self, return_dic=True):
		if(return_dic):
			return self.uni, self.uni_dic
		else:
			return self.uni
	
	def plot_end_states(self, ax=False, n_ats=-1, title="", norm=False, save=False, save_ext='png', save_PATH=''):
		if(n_ats < 0):
			n_ats = len(self.uni_dic)
		to_plot_x = list(self.uni_dic.keys())[:n_ats]
		to_plot_y = list(self.uni[:n_ats,1])
		if(len(self.uni_dic) < n_ats):
			print(n_ats, len(to_plot_x), to_plot_x)
			for i in range(n_ats - len(to_plot_x)):
				to_plot_x.append(' ' * (i + 1))
				to_plot_y.append(0)
			print(to_plot_x)
		if(ax):
			ax.set_title(title)
			if(norm):
				ax.bar(to_plot_x, to_plot_y / self.N)
				ax.set_ylim((0,1))
			else:
				ax.bar(to_plot_x, to_plot_y)
				ax.set_ylim((0,self.N))
			if(self.dim > 8):
				ax.set_xticklabels(to_plot_x, rotation=-45)
			
			return ax
		else:
			plt.figure()
			plt.title(title)
			if(norm):
				plt.bar(list(self.uni_dic.keys())[:n_ats], self.uni[n_ats:,1] / self.N)
			else:
				plt.bar(list(self.uni_dic.keys())[:n_ats], self.uni[:n_ats,1])
			if(self.dim > 8):
				plt.xticks(rotation=90)
			plt.tight_layout()
			
			
			if(save):
				plt.savefig(save_PATH + 'end_states_{:d}cells_initial_{:s}.{:s}'.format(self.N, np.binary_repr(self.s_init, self.dim), save_ext), fmt=save_ext)
		




""" Test run. """
if(__name__ == "__main__"):
	L = [0, 0, 0, 2, 2, 2, 0, 2, 1, -1, -3, -1, 3, 1, -1, 1, 1, 0, -1, 2, -1, -2, -4, -1, 1, -3, -2, 2, 1, -3, -1, 3]
	s_init = '11111'
	genes = ['A', 'B', 'C', 'D', 'E']
	gene_statuses = [1, 1, 1, 1, 1]
	
	N = 1000
	
	sim = Simulation(N, s_init, genes, gene_statuses, L)
	sim.simulate_cells()
	sim.end_states()
	sim.plot_end_states()
	
	plt.show()
	
	
	
	
	
	
