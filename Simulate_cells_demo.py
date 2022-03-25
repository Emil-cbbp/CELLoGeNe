import sys
import os
import CellSimulation as CS
import csv
import copy

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

SMALL_SIZE = 14
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

land = np.array([0, 0, 0, 2, 2, 2, 0, 2, 1, -1, -3, -1, 3, 1, -1, 1, 1, 0, -1, 2, -1, -2, -4, -1, 1, -3, -2, 2, 1, -3, -1, 3])
print(len(land))
labels_land = ['A', 'B', 'C', 'D', 'E']


# Set gene statuses. If 0, the gene cannot change, if 1, it can.
A = 1
B = 1
C = 1
D = 1
E = 1

n = 5	# Number of genes

gene_statuses = [A, B, C, D, E]


N = 10000	# Number of times to simulate an initial cell state	(10000 is reasonable for this small network)


attractor_colour_fig_6a = {int('10110',2): 'blue', int('01010',2): '#00FF00', int('11101',2): 'red', int('11001',2): 'red'}
attr_colour = {'10110': 'blue', '01010': '#00FF00', '11101': 'red', '11001': 'red'}

s_init = [('{:0>' + str(n) + 's}').format(bin(x)[2:]) for x in range(1 << n)]

beta = 1

attr_dist_all = []
for k in range(3):		# Number of times to repeat the experiment
	end_states = []
	end_states_dic = []
	mat = []
	for i,s in enumerate(s_init):
		print(s)
		sim = CS.Simulation(N, s, labels_land, copy.deepcopy(gene_statuses), land, beta)
		sim.simulate_cells()
		sim.end_states()
		es, esd = sim.get_end_states()
		
		res  = [0 if s not in esd else esd[s] for s in s_init]
		mat.append(res)
		end_states.append(es)
		end_states_dic.append(esd)
	mat = np.array(mat)

	

	attr_dist = []
	for a in attr_colour:
		attr_dist.append([0 if a not in end_states_dic[i] else end_states_dic[i][a] for i, s in enumerate(s_init)])
	attr_dist.append(1 - np.sum(attr_dist, axis=0))
	attr_dist = np.array(attr_dist)
	attr_dist_all.append(attr_dist)
	
attr_dist_all = np.array(attr_dist_all)
attr_dist = np.mean(attr_dist_all, axis=0)
attr_dist_std = np.std(attr_dist_all, axis=0)

at = list(attr_colour.keys())



a_d = np.array([attr_dist[0,:], attr_dist[1,:], attr_dist[2,:] + attr_dist[3,:], attr_dist[4,:]])
df = pd.DataFrame(np.transpose(a_d))
df.insert(0, 'states', np.arange(1 << n, dtype=int))

at_ed = at[:2]
at_ed.append(at[2] + '-' + at[3])
at_ed.append('Other')
df = df.rename(columns={x: at_ed[x] for x in range(len(at_ed))})

a_std = np.array([attr_dist_std[0,:], attr_dist_std[1,:], np.sqrt(attr_dist_std[2,:]**2 + attr_dist_std[3,:]**2), attr_dist_std[4,:]])
df_std = pd.DataFrame(np.transpose(a_std))
df_std = df_std.rename(columns={x: at_ed[x] + '_std' for x in range(len(at_ed))})







attr_colour['Other'] = 'grey'
attr_colour[at_ed[2]] = attr_colour[at_ed[2].split('-')[0]]

df_s = df.drop('states', axis=1)
df_q = pd.DataFrame(columns=at_ed)
q = df_s.idxmax(axis=1)
sort = [x for x in range(4)]
for i, a in enumerate(at_ed[:-1]):
	r = np.where(q==a)[0]
	df_q = pd.concat([df_q, df_s.loc[r].sort_values(at_ed[sort[i]])])
if(len(df_q) != 1 << n):
	print('Warning!')
uni = np.unique(q, return_counts=True)
uni_dic = {uni[0][i]: uni[1][i] for i in range(len(uni[0]))}

lab = list(df.columns)

plt.figure(figsize=(10, 6))
x = np.arange(len(df_q))
y = [np.zeros(1<<n)]
u = []
for i, a in enumerate(at_ed):
	y.append(np.array(np.sum(np.array([df_q[at_ed[j]] for j in range(i+1)]), axis=0)))
	plt.fill_between(x, y[i+1], y[i], color=attr_colour[at_ed[i]], label=at_ed[i], alpha=0.9)
for i, a in enumerate(at_ed):
	#plt.fill_between(x, y[i+1] + df_std[a + '_std'], y[i+1], color='black', alpha=0.4)
	#plt.fill_between(x, y[i+1] - df_std[a + '_std'], y[i+1], color='black', alpha=0.4)
	if(a in uni_dic):
		u.append((0 if len(u) == 0 else u[-1])  + uni_dic[at_ed[i]])
u[-1] -= 1
for k, t in enumerate(u):
	plt.vlines(t, -0.04, 1, color='black')
	plt.text((u[k-1] if k > 0 else 0) + (t - (u[k-1] if k > 0 else 0)) / 2, -0.04, at_ed[k], horizontalalignment='center', verticalalignment='top')
		
	

plt.xticks([])
plt.xlabel('States')
plt.ylabel('Probability of reaching attractor')
plt.legend(loc='best')









""" Code for simulating with attractors as initial states. """
N = 10000
attr_colour = {'10110': 'blue', '01010': '#00FF00', '11101': 'red', '11001': 'red'}

beta = 1

s_init = [	'10110', '01010', '11101']


attr_dist_all = []
for k in range(3):
	print('Starting simulation set', k+1)
	end_states = []
	end_states_dic = []
	mat = []
	for i,s in enumerate(s_init):
		print(s)
		sim = CS.Simulation(N, s, labels_land, copy.deepcopy(gene_statuses), land, beta)
		sim.simulate_cells(perturb_init=False)
		sim.end_states()
		es, esd = sim.get_end_states()
		
		res  = [0 if s not in esd else esd[s] for s in s_init]
		mat.append(res)
		end_states.append(es)
		end_states_dic.append(esd)
	mat = np.array(mat)
	
	

	attr_dist = []
	for a in attr_colour:
		attr_dist.append([0 if a not in end_states_dic[i] else end_states_dic[i][a] for i, s in enumerate(s_init)])
	attr_dist.append(1 - np.sum(attr_dist, axis=0))
	attr_dist = np.array(attr_dist)
	attr_dist_all.append(attr_dist)
	
attr_dist_all = np.array(attr_dist_all)
attr_dist = np.mean(attr_dist_all, axis=0)
attr_dist_std = np.std(attr_dist_all, axis=0)

at = list(attr_colour.keys())
a_d = np.array([attr_dist[0,:], attr_dist[1,:], attr_dist[2,:] + attr_dist[3,:], attr_dist[4,:]])
df = pd.DataFrame(np.transpose(a_d))
df.index = s_init

at_ed = at[:2]
at_ed.append(at[2] + '-' + at[3])
at_ed.append('Other')
df = df.rename(columns={x: at_ed[x] for x in range(len(at_ed))})

a_std = np.array([attr_dist_std[0,:], attr_dist_std[1,:], np.sqrt(attr_dist_std[2,:]**2 + attr_dist_std[3,:]**2), attr_dist_std[4,:]])


df_std = pd.DataFrame(np.transpose(a_std))
df_std = df_std.rename(columns={x: at_ed[x] + '_std' for x in range(len(at_ed))})
df_std.index = s_init

keys = list(attr_colour.keys())
attr_colour['Other'] = 'grey'
attr_colour[at_ed[2]] = attr_colour[at_ed[2].split('-')[0]]

map_state_lab = {s_init[i]: at_ed[i] for i in range(len(s_init))}
lab = list(df.columns)

c = [attr_colour[s] for s in s_init]
c.append('grey')

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15,5), sharey=True)
ax = ax.flatten()
for k,s in enumerate(s_init):
	ax[k].bar(np.arange(len(s_init) + 1), np.array(df.loc[s]), yerr=np.array(df_std.loc[s]), color=c)
	ax[k].set_xticks(np.arange(len(lab)))
	ax[k].set_xticklabels(lab, rotation=0)
	ax[k].set_ylim(0,1.05)
	ax[k].grid(axis='y')
	ax[k].set_axisbelow(True)
	
	if('#' in c[k]):
		ax[k].set_title('Initial state: {:s}'.format('Green'))	
	else:
		ax[k].set_title('Initial state: {:s}'.format(c[k].capitalize()))
ax[0].set_ylabel('Probability of reaching attractor')
plt.tight_layout()


