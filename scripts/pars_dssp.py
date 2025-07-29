#!/usr/bin/python
import sys

Norm_Acc={"A" :106.0, "B" :160.0,
			"C" :135.0, "D" :163.0, "E" :194.0,
			"F" :197.0, "G" : 84.0, "H" :184.0,
			"I" :169.0, "K" :205.0, "L" :164.0,
			"M" :188.0, "N" :157.0, "P" :136.0,
			"Q" :198.0, "R" :248.0, "S" :130.0,
			"T" :142.0, "V" :142.0, "W" :227.0,
			"X" :180.0, "Y" :222.0, "Z" :196.0}	

def get_dssp(dsspfile,chain):
	dssp=[]
	s=0
	f=open(dsspfile)
	for line in f:
		if line.startswith('  #  RESIDUE'):
			s=1
			continue
		if s==0: continue
		if line[11] !=chain: continue
		n=int(line[5:10].strip())
		ch=line[11]
		aa=line[13]
		ss=line[16]
		sa=float(line[34:38])
		phi=float(line[103:109])
		psi=float(line[109:115])
		if ss==" ": ss="C"
		dssp.append([n,ch,aa,ss,sa,phi,psi])
	return dssp
			
			
def get_rsa(dssp,msa=Norm_Acc):
	for n,ch,aa,ss,sa,_,_ in dssp:         #for n,ch,aa,ss,sa,phi,psi in dssp
		                                   #in this situation you can use _ to avoid keep them in memory
		if not msa.get(aa,False): continue
		print ("%d\t%s\t%s\t%s\t%.3f" %(n,ch,aa,ss,sa/msa[aa]))
		

	
			

if __name__ =='__main__':
	dsspfile=sys.argv[1]
	chain=sys.argv[2]
	dssp=get_dssp(dsspfile,chain)
	get_rsa(dssp)
	
	
	
	