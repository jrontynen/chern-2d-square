
import os
import re



#For natural sorting:

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]




res = []
for fname in os.listdir("./out"):
	if fname[:5] == 'jobc1':
		if 'SUCCESS' not in open('./out/' + fname).read():
			res.append(fname + " no success")
                #if 'TIME' in open('./out/' + fname).read():
		#	res.append(fname + " time limit")


if len(res) > 0:
	res.sort(key=natural_keys)
	print "\n".join(res)
	





