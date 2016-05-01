from __future__ import division
from itertools import *

class HMM:

	def __init__(self, obs, states, start_p, trans_p, emit_p):

		self.obs = obs
		self.states = states
		self.start_p = start_p
		self.trans_p = trans_p
		self.emit_p = emit_p

	def print_dptable(self, M):

	    print "    ",
	    for i in range(len(M)): print "%12s" % ("%d" % i),
	    print

	    for j in M[0].keys():
	        print "%7s: " %j,
	        for t in range(len(M)):
	            print "%7s" % ("%.10f" % M[t][j]),
	        print

	def forward(self):
		self.alpha = [{}]
		self.prob_alpha = 0

		for i in self.states:
			self.alpha[0][i] = self.start_p[i] * self.emit_p[i][self.obs[0]]
		for t in range(1,len(self.obs)):
			self.alpha.append({})
			for i in self.states:
				self.alpha[t][i] = self.emit_p[i][self.obs[t]]*sum(self.alpha[t-1][j] * self.trans_p[j][i] for j in self.states)
		self.prob_alpha = sum(self.alpha[t][j] for j in self.states)

		return (self.alpha, self.prob_alpha)
		
	def backward(self):
		self.beta = []
		self.prob_beta = 0
		
		for t in range(len(self.obs)):
			self.beta.append({})
			for i in self.states:
				if t == len(self.obs) - 1:
					self.beta[t][i] = 1
				else:
					self.beta[t][i] = 0

		for t in reversed(range(len(self.obs))):
			for i in self.states:
				if t in range(len(self.obs)-1):
					self.beta[t][i] = sum(self.beta[t+1][j] * self.trans_p[i][j] * self.emit_p[j][self.obs[t+1]] for j in self.states)
	
		self.prob_beta = sum(self.start_p[j] * self.emit_p[j][self.obs[0]] * self.beta[0][j] for j in self.states)
		
		return (self.beta, self.prob_beta)
		
	def viterbi(self):
		self.V = [{}]
		self.prob_V = 0
		self.path = {}

		for i in states:
			self.V[0][i] = self.start_p[i] * self.emit_p[i][self.obs[0]]
			self.path[i] = [i]

		for t in range(1,len(self.obs)):
			self.V.append({})
			newpath = {}

			for i in states:
				(self.prob_V, state) = max([(self.V[t-1][j] * self.trans_p[j][i] * self.emit_p[i][self.obs[t]], j) for j in self.states])
				self.V[t][i] = self.prob_V
				newpath[i] = self.path[state] + [i]
	 
			self.path = newpath
	 
		(self.prob_V, state) = max([(self.V[len(self.obs) - 1][i], i) for i in self.states])
		self.path = self.path[state]

		return (self.V, self.path)

	def baum_welch(self,iterations = 1):

			def reestimation():

				self.xi = [] 
				self.gamma = []
				for t in range(len(self.obs)-1):
					self.xi.append({})
					for i in self.states:
						self.xi[t][i] = {}
						for j in self.states:
							first_part = self.alpha[t][i] * self.trans_p[i][j] * self.beta[t+1][j] * self.emit_p[j][self.obs[t+1]]
							second_part = sum(sum(self.alpha[t][i] * self.trans_p[i][j] * self.beta[t+1][j]\
							 * self.emit_p[j][self.obs[t+1]] for j in self.states) for i in self.states)
							self.xi[t][i][j] = first_part / second_part
							
				for t in range(len(self.obs)-1):
					self.gamma.append({})
					for i in states:
						first_part = self.alpha[t][i] * self.beta[t][i]
						second_part = sum(self.alpha[t][j] * self.beta[t][j] for j in states)
						self.gamma[t][i] = first_part / second_part

				start_p_updated = {}
				for i in states:
					start_p_updated[i] = round(self.gamma[0][i],2)

				trans_p_updated = {}
				for i in states:
					trans_p_updated[i] = {}
					for j in states:
						first_part = sum(self.xi[t][i][j] for t in range(len(self.obs)-1))
						second_part = sum(self.gamma[t][i] for t in range(len(self.obs)-1))
						trans_p_updated[i][j] = round(first_part / second_part,2)

				emit_p_updated = {}
				for j in self.states:
					emit_p_updated[j] = {}
					for k in self.obs:
						first_part = 0
						second_part = 0
						for t in range(len(self.obs)-1):
							if k == self.obs[t]: 
								first_part += self.gamma[t][j]
							second_part += self.gamma[t][j]
						emit_p_updated[j][k] = round(first_part / second_part,2)
						

				self.start_p = start_p_updated
				self.trans_p = trans_p_updated
				self.emit_p = emit_p_updated

			for i in range(iterations):
				HMM.forward()
				HMM.backward()
				reestimation()

	def execute(self,*args):
		for arg in args:
			if arg == 'forward':
				print 'Forward algoritms'
				self.forward()
				self.print_dptable(self.alpha)
				print self.prob_alpha, observations
				print('\n')
			if arg == 'backward':
				print 'Backward algoritms'
				self.backward()
				self.print_dptable(self.beta)
				print self.prob_beta, observations
				print('\n')
			if arg == 'viterbi':
				print 'Viterbi algoritms'
				self.viterbi()
				self.print_dptable(self.V)
				print self.prob_V, self.path, observations
				print('\n')
			if arg == 'baum_welch':
				print 'Baum-Welch algoritms'
				self.forward()
				self.backward()	
				self.baum_welch()
				print 'transition_probability =',self.trans_p
				print 'emission_probability =',self.emit_p
				print 'start_probability =',self.start_p
				print('\n')

states = ('Healthy', 'Fever')

observations = ('Normal', 'Cold', 'Dizzy', 'Normal')

start_probability = {'Healthy' : 0.6, 'Fever' : 0.4}

transition_probability = {
    'Healthy' : {'Healthy' : 0.7, 'Fever' : 0.3},
    'Fever' : {'Healthy' : 0.4, 'Fever' : 0.6},
}

emission_probability = {
    'Healthy' : {'Normal' : 0.5, 'Cold' : 0.4, 'Dizzy' : 0.1},
    'Fever' : {'Normal': 0.1, 'Cold' : 0.3, 'Dizzy' : 0.6}
}
		
HMM = HMM(observations, states, start_probability, transition_probability, emission_probability)

HMM.execute('forward')
HMM.execute('backward')
HMM.execute('viterbi')
HMM.execute('baum_welch')
HMM.execute('forward')
