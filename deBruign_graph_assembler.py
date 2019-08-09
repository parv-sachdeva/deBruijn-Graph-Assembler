import numpy as np
import pandas as pd
import re

class DNA():
	def __init__(self, seq):
		self.seq=seq
		self.seq_len=len(seq)

	def all_kmer(self, k):
		self.kmers=[]
		for i in range(self.seq_len):
			if i+k <= self.seq_len:
				self.kmers.append(self.seq[i:k+i])
		return self.kmers

	def all_k1mer(self, kmers=None):
		if not kmers:
			tmp=self.kmers
		else:
			tmp=kmers
		self.k1mers=[]
		for kmer in tmp:
			self.k1mers.append(kmer[:-1])
			self.k1mers.append(kmer[1:])
		return self.k1mers


	def unique_kmers(self, kmers=None):
		'''Find unique kmers, to find kmers run all_kmer(k) first!'''
		if kmers:
			return set(kmers)
		else:
			return set(self.kmers)

class ShortestCommonSuperstring(DNA):
	def __init__(self, kmers=None):
		if kmers!=None:
			self.kmers=kmers
			self.scs()

		else:
			print('Warning! No kmers provided. You can load sequences using load_seq() function.')

	def load_seq(self, seq, k):
		self.k=k
		self.seq=seq
		super().__init__(seq)
		self.kmers=super().all_kmer(k)
		self.k1mers=super().all_k1mer()
		#self.scs()

	def overlap(self, str1, str2):
		'''Finding overlap in suffix of str1 and prefix of str1'''
		len1=len(str1)
		len2=len(str2)
		combined=""
		for i in range(len1):
			if re.match(str1[i:], str2):
				offset=re.match(str1[i:], str2).span()[1]
				combined=str1+str2[offset:]
				return([offset, combined])
				break;

	def graph(self):

		#self.uniq_kmers=super().unique_kmers(kmers=self.kmers) # vertices
		self.edges=[]	#
		tmp=[]
		for kmer1 in self.kmers:
			for kmer2 in self.kmers:
				if kmer1!=kmer2:
					tmp=self.overlap(kmer1, kmer2)
					if not tmp:
						continue
					elif tmp[0]==0:
						continue
					else:
						self.edges.append((kmer1, kmer2, tmp[0], tmp[1]))	##combined=str1+str2[offset:]

		#print(self.uniq_kmers, self.edges)

	def scs(self):
		'''Recursive function for creating shortest common superstring'''
		self.graph()
		self.max=0
		for edge in self.edges:
			#print(edge)
			if self.max< edge[2]:
				self.max=edge[2]
		i=len(self.kmers);tmp=[]
		for edge in self.edges:
			if edge[2]==self.max:
				tmp=edge
				#del self.edges[i]
				#print(self.kmers)

				self.kmers.remove(edge[0])
				self.kmers.remove(edge[1])
				self.kmers.append(tmp[3])
				#print(self.kmers)
				self.scs()
				break
		if len(self.kmers)>=i:
			print('Shortest common superstring : ', self.kmers)
			return(self.kmers)
		i=len(self.kmers)
					
class deBruijnGraph(ShortestCommonSuperstring):
		def __init__(self, kmers=None):
			if kmers!=None:
				self.kmers=kmers
				self.k1mers=super().all_k1mer(self.kmers)

			else:
				print('Warning! No kmers provided. You can load sequences using load_seq() function.')

		def dgraph(self):
			'''Function to create a deBruign graph'''
			i=0
			self.glist=[[self.k1mers[i], self.k1mers[i+1]] for i in range(0, len(self.k1mers), 2)]
			self.edges=[]
			for pair1 in self.glist:
				for pair2 in self.glist:
					if pair1==pair2:
						continue
					else:
						#print(pair1, pair2)
						tmp=super().overlap(pair1[1], pair2[0])
						#print(tmp)
						if not tmp:
							continue
						elif tmp[0]==0:
							continue
						else:
							self.edges.append((pair1, pair2, tmp[0], tmp[1]))
			#print(self.edges)

		def assemble(self, edges=None, glist=None):
			'''Recursive function to traverse through deBruign graph'''
			self.dgraph()
			self.assembly=""
			if edges!=None:
				self.edges=edges
			if glist!=None:
				self.glist=glist
			max=0
			i=len(self.kmers)
			for edge in self.edges:
				if edge[2] > max:
					max=edge[2]
					#print(edge[0][0])
			for edge in self.edges:
				if edge[2]==max:
					self.assembly+=edge[0][0][0]+edge[0][1]+edge[1][0][edge[2]:]+edge[1][1][-1]
					self.edges.remove(edge)
					#print(self.assembly)
					self.kmers.remove(edge[0][0][0]+edge[0][1])
					self.kmers.remove(edge[1][0][0]+edge[1][1])
					self.kmers.append(self.assembly)
					self.k1mers=super().all_k1mer(self.kmers)
					#update kmers and k1mers
					self.assemble()
					break
				#print(self.kmers)
			if len(self.kmers)>=i:
				print('Assembled string : ', self.kmers)
				return(self.kmers)
			i=len(self.kmers)

