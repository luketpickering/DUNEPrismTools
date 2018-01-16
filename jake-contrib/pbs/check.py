import os
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--setting', type=str, help='which dir', default='FHC')
parser.add_argument('--subdir', type=str, help='which subdir', default='')
print parser.parse_args().setting

directory = 'sub_lists/' 
setting = parser.parse_args().setting

#setting = 'FHC_FD'

new_sub_files = []
files = 0
extra = 0
newsub = open(setting + '_new_submit_0.list','w')
subdir = setting + parser.parse_args().subdir
#output = listdir('/home/calcuttj/DUNEPrismSim/' + subdir)
output = [] 
for f in os.listdir('/home/calcuttj/DUNEPrismSim/' + subdir):
  p = os.path.join('/home/calcuttj/DUNEPrismSim/' + subdir,f)
  if os.path.isdir(p):
    continue
  output.append(f)

#submit = listdir(subdir)
submit = []
for f in os.listdir(subdir):
  p = os.path.join(subdir,f)
  if os.path.isdir(p):
    continue
  submit.append(f)

strip_output = []
strip_submit = []

for f in output:
  #print f.split('.')
  strip_output.append('.'.join((f.split('.')[0:2])))
  print '.'.join((f.split('.')[0:2]))
  #print f

good = 0
bad = 0
pre = '.'.join(submit[0].split('.')[0:2])

for f in submit:
  #print f.split('.')
  strip_submit.append('.'.join(f.split('.')[2:4]))
  print '.'.join(f.split('.')[2:4])
  if '.'.join(f.split('.')[2:4]) not in strip_output:
#    print int(f.split('.')[3]), "missing" 
    files = files + 1 
    bad = bad + 1 
    if (files > 250):
      newsub.close()
      extra = extra + 1 
      newsub = open(setting + '_new_submit_' + str(extra) + '.list','w')
      files = 1
    newsub.write(subdir + pre + '.' + '.'.join(f.split('.')[2:4]) + '.pbs\n')
  else:
    good = good + 1 
 
print good, bad 
#print pre



