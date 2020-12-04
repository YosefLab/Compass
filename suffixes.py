def suffix_array(S):
	indices = list(range(1,len(S)+1))
	return sorted(indices, key=lambda i:S[i-1:])
	
def suffix_array_zind(S):
	indices = list(range(len(S)))
	return sorted(indices, key=lambda i:S[i:])
	
def rotations(S):
	suffixes = []
	for i in range(len(S)):
		suffixes += [S[i:]+S[:i]]
	return sorted(suffixes)

def BWT(S):
	return [s[-1] for s in rotations(S)]