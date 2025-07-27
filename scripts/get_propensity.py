#!/usr/bin/python
import sys
import numpy as np

def get_seq(seqfile):
	d={}
	f=open(seqfile)
	for line in f:
		if line.startswith('>'):
			pid=line.strip()[1:].split(':')
			sid=pid[0]+':'+pid[1]
			d[sid]=''
		else:
			d[sid]=d[sid]+line.strip()
	return d
	


def get_propensity(dseq,dss,s='H',aa='ACDEFIKLMNPQRSTVWY'):
	naa=np.zeros(20)         #number of the aminoacids
	nas=np.zeros(20)			#number of the aminoacid structures
	ps=np.zeros(20)
	ks=list(dseq.keys())    #keys
	for k in ks:
		seq=dseq.get(k,'')
		ss=dss.get(k,'')	
		if len(ss)!=len(seq) or len(seq)==0: continue
		for i in range(len(seq)):
			p=aa.find(seq[i])
			if p>-1: naa[p]+=1
			if ss[i]==s: nas[p]+=1
		n=float(naa.sum())
		ns=float(nas.sum())
		pss=ns/n
		paa=naa/n
		pas=nas/n
		for i in range(len(aa)):
			ps[i]=pas[i]/(pss*paa[i])
			print (f"{aa[i]}\t{ps[i]}")
		return ps


if __name__ == '__main__':
	seqfile=sys.argv[1]
	ssfile=sys.argv[2]
	s=sys.argv[3]
	dseq=get_seq(seqfile)
	dss=get_seq(ssfile)
	get_propensity(dseq,dss,s)
	
	
	
		