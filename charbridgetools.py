#!/usr/bin/python
import Levenshtein

def bridgehunter(bridge,line):
	position,bridgenum,orientation=bridgepos_multicheck(bridge,line)

	if bridgenum == 0:
		orientation = 'F'
		rna = ''
		dna = ''
		position = -99
		NNN = ''  # VCM added
		nnngate = 0  # VCM added
		dnagate = 0  # VCM added
	elif bridgenum == 1:
		#jackpot!
		#is it fwd or reverse?
		if orientation=='F':
			#forward
			prna=line[0:position]
			pdna=line[(position+len(bridge)):]

			#rna,dna=trimrnadna(prna,pdna)
			rna,dna,NNN,nnngate,dnagate=trimrnadna(prna,pdna)  # VCM added
		elif orientation=='R':
			#reverse
			#so we want to reverse the sequences
			backwardsdna=line[0:position]
			backwardsrna=line[(position+len(bridge)):]

			pdna=reverse_complement(backwardsdna)
			prna=reverse_complement(backwardsrna)

			#rna,dna=trimrnadna(prna,pdna)
			rna,dna,NNN,nnngate,dnagate=trimrnadna(prna,pdna)  # VCM added
			position=len(line)-position-len(bridge)-1 #so it equals the equivalent if fwd
		else:
			print "no orientation information provided"
	elif bridgenum >= 1:
		orientation='F'
		rna=''
		dna=''
		position=-100
		NNN=''  # VCM added
		nnngate=0  # VCM added
		dnagate=0  # VCM added
	else:
		pass
	# return rna,dna,orientation,position,bridgenum
	return rna,dna,orientation,position,bridgenum,NNN,nnngate,dnagate  # VCM added

def trimrnadna(prna,pdna):
	##### Added code below #####
	# assumed: NNN = 'AANNN' + core bridge == 'AAACCGGCGTCCAAG' + dpnii 'GATC'
	NNN = prna[-5:]  # prna[-3:] == 'CAT' from 'AACAT'
	cleanrna = prna[:-5]  # remove rna_end_bridge = 'AANNN'
	nnngate = 1
	#if NNN[:2] == 'AA':  ## the two 3' AA's are included in the core bridge sequence in char_bridge_trackall
	#	nnngate = 1  # if random trimer NNN is flanked by AA on 5' side
	#else:
	#	nnngate = 0  # if random trimer NNN is NOT flanked by AA on both sides then NNN might not be the random trimer.
	##### Added code above #####

	cleandna=pdna
	##### Added code below #####
	bridge_tail_b4_ligation = 'GATCTTTAATTAAGTCGCAG'  # 'GATCTTTAATTAAGTCGCAGATC'  # for WT with Pac1
	#bridge_tail_b4_ligation = 'GATCTGCGGCCGCGTCGCAG'  # 'GATCTGCGGCCGCGTCGCAGATC'  # for KD with Not1
	dpnii = 'GATC'
	if fuzz_align_once(bridge_tail_b4_ligation, cleandna, 1):
		#cleandna="GATC"+pdna[len(bridge_tail_b4_ligation):]  # VCM edit  # position = 0
		# cleandna="GATC"+pdna[(len(i)-1):]  # default  # VCM: i'm pretty sure the -1 does not remove the C in GATC
		# cleandna=pdna[(len(i)):]  # default
		dnagate = 0
	elif dpnii == cleandna[:len(dpnii)]:
		dnagate = 1
	else:
		# if there is no gatc at the 5' end of the DNA do not allow it to be aligned
		dnagate = 0
	##### Added code above #####

	return cleanrna,cleandna,NNN,nnngate,dnagate


def fuzz_align_once(s_seq,l_seq,mismatch):
	for i, base in enumerate(l_seq):  # loop through equal size windows
		l_subset = l_seq[i:i+len(s_seq)]
		dist = Levenshtein.distance(l_subset, s_seq)
		if dist <= mismatch:  # find first then break
			return i, dist
		#l_subset = l_seq[i:i+len(s_seq)]
		#if mismatch == 0:  #  VCM added
		#	if l_subset == s_seq:  # find first then break  #  VCM added
		#		return i, 0  #  VCM added
		#else:
		#	dist = Levenshtein.distance(l_subset, s_seq)
		#	if dist <= mismatch:  # find first then break
		#		return i, dist


def bridgepos_multicheck(bridge, line):
	count=0
	#try finding the bridge
	fuzz=fuzz_align_once(bridge,line,1)
	#if you didnt find it, try the revcomp of the bridge
	if not(fuzz):
		fuzz_rc=fuzz_align_once(reverse_complement(bridge),line,1)
		if not(fuzz_rc):
			#if no bridge either way, just give up.
			return 0,0,'F'
		else:
			#if we found the bridge backwards once, check to make sure we don't see it again
			restofline=line[(fuzz_rc[0]+len(bridge)):]
			fuzz_rc_2=fuzz_align_once(reverse_complement(bridge),restofline,1)
			if not(fuzz_rc_2):
				#so this means we found it once but not twice
				return fuzz_rc[0],1,'R'
			else:
				#this means we found it multiple times
				return fuzz_rc[0],2,'R'
	else:
		#so we found it in the fwd direction
		#check the rest of the read in the fwd direction
		restofline=line[(fuzz[0]+len(bridge)):]
		fuzz_fw=fuzz_align_once(bridge,restofline,1)
		#and the rc
		fuzz_rc=fuzz_align_once(reverse_complement(bridge),restofline,1)
		#if either returns something, it means there's a double bridge
		if fuzz_fw or fuzz_rc:
			return fuzz[0],2,'F'
		else:
			#if there isn't, return what we want
			return fuzz[0],1,'F'

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):	
	for k,v in alt_map.iteritems():
		seq = seq.replace(k,v)
	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	for k,v in alt_map.iteritems():
		bases = bases.replace(v,k)
	return bases

