import os, sys
try:
  from iotbx import pdb
  from libtbx import easy_run
  from __future__ import division
except ImportError:
  sys.stderr.write("This program requires the Phenix environment to run.\n")
  sys.stderr.write("Please source the Phenix environment (phenix/build/setpaths.sh) and try again.\n")
  sys.exit()

#todo list
#zstrack weighted balls
#topout at ~0.6, bottomout large enough to draw near-misses underneath
#consider weighted spheres for vdw and clash
#subprocess instead of os.system

#--- Get file from user, other setup
try:
  pdbfile = sys.argv[1]
except IndexError:
  sys.stderr.write("""
Please supply a single PDB-formatted file.
Add Hydrogens to your file with Reduce if it does not already contain H's.
mmCIF-formatted files may be converted to PDB with phenix.pdb_as_cif.
  (Do this before running Reduce)

Output will be a kinemage printed to stdout.
This kinemage is a sequence contact map with categories from probe contact types
  wc = vdw wide contact
  cc = vdw close contact
  hb = hydrogen bond
  so = small overlap
  bo = bad overlap (clash)
If the input file has multiple chains, chains are separated by white lines in
  the output.

""")
  sys.exit()

def make_one_view(name="Unnamed view", span="100", center="0 0 0", view_num="1"):
  sys.stdout.write("@"+str(view_num)+"viewid {"+name+"}\n")
  sys.stdout.write("@"+str(view_num)+"span "+str(span)+"\n")
  sys.stdout.write("@"+str(view_num)+"zslab 200.0\n")
  sys.stdout.write("@"+str(view_num)+"center "+center+"\n")
  sys.stdout.write("@"+str(view_num)+"matrix 1 0 0 0 1 0 0 0 1\n")

def make_all_views(chain_indices, chain_ids):
  end = chain_indices[-1]
  view_num = 1
  base_span = end
  span = base_span + base_span/20
  xcenter = end/2
  ycenter = end/2
  center = str(xcenter)+" "+str(ycenter)+" 0"
  make_one_view(name="overview", span=span, center=center, view_num=view_num)
  start_index = 0
  end_index = 1
  while start_index < len(chain_indices)-1:
    view_num = end_index+1
    base_span = chain_indices[end_index] - chain_indices[start_index]
    span = base_span + base_span/20
    xcenter = (chain_indices[end_index] + chain_indices[start_index])/2
    ycenter = (chain_indices[end_index] + chain_indices[start_index])/2
    center = center = str(xcenter)+" "+str(ycenter)+" 0"
    make_one_view(name="chain "+chain_ids[start_index], span=span, center=center, view_num=view_num)
    start_index += 1
    end_index += 1

#@1viewid {Unnamed view}
#@1span 868.32715
#@1zslab 200.0
#@1center 307 307 0
#@1matrix 1 0 0 0 1 0 0 0 1
#@2viewid {this is a view}
#@2span 349.60812
#@2zslab 200.0
#@2center 419 422 0
#@2matrix 1 0 0 0 1 0 0 0 1


#--- Load into phenix, set up numbering conversions
pdb_io = pdb.input(pdbfile)
hierarchy = pdb_io.construct_hierarchy()

indexing = {}
i = 2
all_indices = []
chain_indices = [0]
chain_ids = []

#Each residue gets a unique index that becomes its x/y coordinate
#Some extra space is left between chains for visual clarity
for chain in hierarchy.chains():
  print chain.id, chain.is_protein(), chain.is_na()
  if not (chain.is_protein() or chain.is_na()): continue
  for rg in chain.residue_groups():
    resid = chain.id.strip()+rg.resseq
    indexing[resid] = i
    all_indices.append(i)
    i+=1
  i+=1 # gap between chains
  chain_indices.append(i)
  chain_ids.append(chain.id)
  i+=2


#probe_out = easy_run.fully_buffered(probe_command, stdin_lines=input_str)
#return probe_out.stdout_lines

#constants for contact type hierarchy
#lower values have precedence
HB = 0
BO = 1
SO = 2
CC = 3
WC = 4

#---Run Probe
probecommand = 'phenix.probe -quiet -u -condensed -self -mc -nohets -nowaters ALL '+pdbfile
probe_out = easy_run.fully_buffered(probecommand)
probefile = probe_out.stdout_lines
#os.system(probecommand)

#probefile = open("temp_probe_out.probe")

#contact type to color and size mappings are stored here
#technically, these are balls not dots.  Oh well.
dottypes = ['wc','cc','hb','hb_weighted','so','bo']
dotcolor = {
  'wc':'blue',
  'cc':'green',
  'hb':'sea',
  'hb_weighted':'sea',
  'so':'yellow',
  'bo':'hotpink'}
dotsize = {
  'wc':'0.5',
  'cc':'0.45',
  'hb':'0.4',
  'hb_weighted':'1',
  'so':'0.35',
  'bo':'0.3'}

mcatoms = [' N  ',' CA ',' C  ',' O  ',' H  ']
#DNA
mcatoms += [' P  ',' OP1',' OP2'," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"]
mcatoms += [" H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"]
#RNA
mcatoms += [" O2'","HO2'"]

kinglines = {}
for dottype in dottypes:
  kinglines[dottype+'_mcmc'] = []
  kinglines[dottype+'_mcsc'] = []
  kinglines[dottype+'_scsc'] = []

#----------Begin condensed flag management----------
#The onedoteach flag has some issues with missing h-bonds that are a little too
#  close.  Those bonds get their one dot assigned as bo or so, not hb
#This version of the script uses the -condensed flag, but must further condense
#  the probe output into one dot per contact
onedoteach = {}
for line in probefile:
  #:1->1:so: A 513 CYS  SG  : A 509 VAL  O   :0.00:-0.075:-0.075:-0.176:39.885:60.076:0.038:-0.0236:S:O:-0.152:39.914:60.078:44.60:36.72
  x = line.split(":")
  dottype = x[2]
  if dottype == 'hb': dottype = HB
  elif dottype == 'bo': dottype = BO
  elif dottype == 'so': dottype = SO
  elif dottype == 'cc': dottype = CC
  else: dottype = WC
  src = x[3]
  trg = x[4]
  contactid = src+trg
  if contactid not in onedoteach:
    onedoteach[contactid] = {'dottype':dottype,'line':line}
  elif onedoteach[contactid]['dottype'] > dottype:
    onedoteach[contactid] = {'dottype':dottype,'line':line}

for item in onedoteach.values():
  line = item['line']
#----------End condensed flag management----------
#for line in probefile:
  #:1->1:so: A 513 CYS  SG  : A 509 VAL  O   :0.00:-0.075:-0.075:-0.176:39.885:60.076:0.038:-0.0236:S:O:-0.152:39.914:60.078:44.60:36.72 ode
  #:1->1:wc: A 192 GLY  CA  : A 332 ILE  O   :10  :0.199 :0.383 :32.597:61.654:15.605:0.000:0.0060 :C:O:32.597:61.654:15.605:18.64:15.92 con
  x = line.split(":")
  dottype = x[2]
  src = x[3]
  trg = x[4]
  srcchain = src[0:2].strip()
  srcnum = src[2:6]
  srcatom = src[11:15]
  trgchain = trg[0:2].strip()
  trgnum = trg[2:6]
  trgatom = trg[11:15]
  srcid = srcchain+srcnum
  trgid = trgchain+trgnum
  srcindex = str(indexing[srcid])
  trgindex = str(indexing[trgid])
  #assemble on-click text for the dot
  kingtext = ":".join([dottype,src,trg])
  #assemble dotlist line, using src and trg indices are x and y coordinates
  pairtype = ''
  if srcatom in mcatoms and trgatom in mcatoms:
    pairtype = '_mcmc'
  elif srcatom not in mcatoms and trgatom not in mcatoms:
    pairtype = '_scsc'
  else:
    pairtype = '_mcsc'
  if dottype == 'hb':
    mingap = float(x[6]) #this works for both con and ode
    ball_radius = mingap/-0.6/2
    if ball_radius > 0.5:
      ball_radius = 0.5
    elif ball_radius < 0.1:
      ball_radius = 0.1
    z = (0.5-ball_radius)
    #kinglines['hb_weighted'+pairtype].append("{"+kingtext+" : "+str(mingap)+"} r="+str(ball_radius)+" "+srcindex+" "+trgindex+" "+str(z)+" r="+str(ball_radius))
    kinglines["hb_weighted"+pairtype].append("{"+kingtext+" : "+str(mingap)+"} r=%.3f %s %s %.3f"  % (ball_radius, srcindex, trgindex ,z))
  kinglines[dottype+pairtype].append("{"+kingtext+"} "+srcindex+" "+trgindex)

sys.stdout.write("@kinemage {"+os.path.basename(pdbfile)+"}\n")
sys.stdout.write("@flat\n") #start in flatland
sys.stdout.write("@hsvcolor {lightgray} 0 0 35 0 0 80\n")
#This is a custom color defintition for king,
#  formatted as a hue saturation value triplet for black backgound,
#  then an optional second hsv triplet for white background

make_all_views(chain_indices, chain_ids)

#print the frame
sys.stdout.write("@group {chains}\n")
#the light gray lines enclose spaces denoting residue pairs
sys.stdout.write("@vectorlist {residues} color=lightgray off\n")
end = str(chain_indices[-1])
for resindex in all_indices:
  sys.stdout.write("{} P 0 "+str(resindex-0.5)+"\n")
  sys.stdout.write("{} "+end+" "+str(resindex-0.5)+"\n")
  sys.stdout.write("{} P "+str(resindex-0.5)+" 0\n")
  sys.stdout.write("{} "+str(resindex-0.5)+" "+end+"\n")
#there's a little extra step to sys.stdout.write(the last line for each chain
for chainindex in chain_indices:
  if chainindex == 0: continue
  sys.stdout.write("{} P 0 "+str(chainindex-1.5)+"\n")
  sys.stdout.write("{} "+end+" "+str(chainindex-1.5)+"\n")
  sys.stdout.write("{} P "+str(chainindex-1.5)+" 0\n")
  sys.stdout.write("{} "+str(chainindex-1.5)+" "+end+"\n")
#the white (default) lines separate each chain and enclose the whole
sys.stdout.write("@vectorlist {borders} width=2\n")
for chainindex in chain_indices:
  sys.stdout.write("{} P 0 "+str(chainindex)+"\n")
  sys.stdout.write("{} "+end+" "+str(chainindex)+"\n")
  sys.stdout.write("{} P "+str(chainindex)+" 0\n")
  sys.stdout.write("{} "+str(chainindex)+" "+end+"\n")


#print the actual content
#print order matters so that balls stack on each other from large to small
for dottype in dottypes:
  sys.stdout.write("@group {"+dottype+"} animate\n")
  sys.stdout.write("@subgroup {Mc Mc} dominant master= {all Mc Mc}\n")
  sys.stdout.write("@balllist {"+dottype+"_mcmc} color="+dotcolor[dottype]+" nohighlight radius="+dotsize[dottype]+"\n")
  for dot in kinglines[dottype+"_mcmc"]:
    sys.stdout.write(dot+"\n")
  sys.stdout.write("@subgroup {Mc Sc} dominant master= {all Mc Sc}\n")
  sys.stdout.write("@balllist {"+dottype+"_mcsc} color="+dotcolor[dottype]+" nohighlight radius="+dotsize[dottype]+"\n")
  for dot in kinglines[dottype+"_mcsc"]:
    sys.stdout.write(dot+"\n")
  sys.stdout.write("@subgroup {Sc Sc} dominant master= {all Sc Sc}\n")
  sys.stdout.write("@balllist {"+dottype+"_scsc} color="+dotcolor[dottype]+" nohighlight radius="+dotsize[dottype]+"\n")
  for dot in kinglines[dottype+"_scsc"]:
    sys.stdout.write(dot+"\n")

