import os, sys
try:
  from iotbx import pdb
  from libtbx import easy_run
  #from iotbx.pdb import secondary_structure
except ImportError:
  sys.stderr.write("This program requires the Phenix environment to run.\n")
  sys.stderr.write("Please source the Phenix environment (phenix/build/setpaths.sh) and try again.\n")
  sys.exit()
#parse probe output for dotcoutns and # of hbonds of different strengths

#--- Get file from user, other setup
try:
  pdbfile = sys.argv[1]
except IndexError:
  sys.stderr.write("""
Please supply a single PDB-formatted file.
Add Hydrogens to your file with Reduce if it does not already contain H's.
mmCIF-formatted files may be converted to PDB with phenix.pdb_as_cif.
  (Do this before running Reduce)
""")
  sys.exit()

def find_resolution(pdbfilepath):
  pdbfile = open(pdbfilepath)
  for line in pdbfile:
    if line.startswith('ATOM'):
      return None
    elif line.startswith('REMARK   2 RESOLUTION'):
      return line.split()[-2]
  return None

#mainchain atom names
#Protein
mcatoms = [' N  ',' CA ',' C  ',' O  ',' H  ']
#DNA
mcatoms += [' P  ',' OP1',' OP2'," O5'"," C5'"," C4'"," O4'"," C3'"," O3'"," C2'"," C1'"]
mcatoms += [" H5'","H5''"," H4'"," H3'"," H2'","H2''"," H1'"]
#RNA
mcatoms += [" O2'","HO2'"]

try:
  from iotbx import pdb
except ImportError:
  phenix_env = "probe"

pdbfile = sys.argv[1]

probe_command = "phenix.probe -u -condense -self -mc -quiet -NOVDWOUT -NOCLASHOUT ALL "+pdbfile


probe_out = easy_run.fully_buffered(probe_command)
probe_results = probe_out.stdout_lines

#probe_file = open("probe_parsing_temp_file.probe")

counts = {'mcmc':0,'mcsc':0,'scsc':0,'het':0,'water':0}
mmpair_counts = {
'mcmc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'mcsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'scsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'het':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'water':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0}
}

dotcounts = {'mcmc':0,'mcsc':0,'scsc':0,'het':0,'water':0}
mmpair_dotcounts = {
'mcmc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'mcsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'scsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'het':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'water':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0}
}

strongcounts = {'mcmc':0,'mcsc':0,'scsc':0,'het':0,'water':0}
mmpair_strongcounts = {
'mcmc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'mcsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'scsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'het':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'water':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0}
}

midcounts = {'mcmc':0,'mcsc':0,'scsc':0,'het':0,'water':0}
mmpair_midcounts = {
'mcmc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'mcsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'scsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'het':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'water':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0}
}

weakcounts = {'mcmc':0,'mcsc':0,'scsc':0,'het':0,'water':0}
mmpair_weakcounts = {
'mcmc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'mcsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'scsc':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'het':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0},
'water':{'protein':0, 'na':0, 'protein-na':0, 'protein-het':0, 'na-het':0, 'protein-water':0, 'na-water':0}
}

PROTEIN = 0
NA = 1
HET = 2
WATER = 3

residue_counts = [0,0,0,0] #keyed by indices above

#--- Load into phenix, set up numbering conversions
pdb_io = pdb.input(pdbfile)
hierarchy = pdb_io.construct_hierarchy()
moltypes = {}
for chain in hierarchy.chains():
  mmtype = 0
  if chain.is_protein(): moltype = PROTEIN
  elif chain.is_na(): moltype = NA
  else: moltype = HET
  for rg in chain.residue_groups():
    for ag in rg.atom_groups():
      if ag.resname == "HOH": moltype = WATER
      #probe line: A  67 TRP  CH2 : A  75 TRP  HE1 :
      #build resid as ccnnnnirrra: | A  67 TRP |
      resid = "%2s%4s%1s%3s%1s" % (chain.id,rg.resseq,rg.icode,ag.resname,ag.altloc)
      moltypes[resid] = moltype
      residue_counts[moltype]+=1

for line in probe_results:
    #Probe Unformatted Output:
    #name:pat:type:srcAtom:targAtom:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    #for condensed output we have:
    #name:pat:type:srcAtom:targAtom:*dotcount*:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    ###'name' is set by the user on the command line
    ###'pat' is one of 1->1, 1->2, or 2->1; where 1 is src and 2 is targ.
    ###'type' is one of wc, cc, so, bo, hb (wide/close contact, small/bad overlap, h-bond).
    ###'srcAtom' and 'targAtom' follow the pattern CNNNNITTT AAAAL, where C is chain, N is number, I is insertion code, T is residue type, A is atom name, and L is alternate conformation flag.
    ###'*dotcount*' is condensed-output-only, and gives the number of dots in the contact
    ###'min-gap' is the distance between atoms, minus their van der Waals radii; i.e., the distance of closest approach for their vdW surfaces. gap is the distance between vdW surfaces at the current dot. Negative values indicate overlap (clashes or H-bonds).
    ###'x','y','z' is a point on the vdW surface; 'spX','spY','spZ' is tip of spike, if any (same as x,y,z for contacts)
    ###'score' is "this dot's contribution to the [Probe] score" (scaled already? YES)
    ###'stype' and 'ttype' are heavy-atom element name (C, N, O, etc)

  if not line.strip(): continue #averts an IndexError problem with empty lines
  x = line.split(':')

  name = x[0]
  pattern = x[1]
  interactiontype = x[2]
  if not interactiontype == 'hb': continue #skip non-h-bonds

  src = x[3]
  srcChain =    src[0:2].strip()
  srcNum =      int(src[2:6].strip())
  srcIns =      src[6:7]#.strip()
  srcResname =  src[7:10].strip()
  if srcResname == 'HOH': continue #skip waters
  srcAtom = src[11:15]#.strip()
  srcAlt =      src[15:16].strip()
  srcResid = src[0:10]+src[15:16]

  trg = x[4]
  #going to count dots per bond as a measure of strength instead
  trgChain =    trg[0:2].strip()
  trgNum =      int(trg[2:6].strip())
  trgNumStr =   trg[2:6]
  trgIns =      trg[6:7]#.strip()
  trgResname =  trg[7:10].strip()
  trgAtom = trg[11:15]#.strip()
  trgAlt =      trg[15:16].strip()
  trgResid = trg[0:10]+trg[15:16]

  dotcount = x[5]
  mingap = float(x[6])

  src_moltype = moltypes[srcResid]
  trg_moltype = moltypes[trgResid]

  pairtype = ''
  if srcAtom in mcatoms and trgAtom in mcatoms:
    pairtype = 'mcmc'
  elif srcAtom not in mcatoms and trgAtom not in mcatoms:
    pairtype = 'scsc'
  else:
    pairtype = 'mcsc'
  if src_moltype == WATER or trg_moltype == WATER:
    pairtype = 'water'
  elif src_moltype == HET or trg_moltype == HET:
    pairtype = 'het'

  if moltypes[srcResid] == PROTEIN:
    if moltypes[trgResid] == PROTEIN:
      mmpair = 'protein'
    elif moltypes[trgResid] == NA:
      mmpair = 'protein-na'
    elif moltypes[trgResid] == HET:
      mmpair = 'protein-het'
    else:
      mmpair = 'protein-water'
  elif moltypes[srcResid] == NA:
    if moltypes[trgResid] == PROTEIN:
      mmpair = 'protein-na'
    elif moltypes[trgResid] == NA:
      mmpair = 'na'
    elif moltypes[trgResid] == HET:
      mmpair = 'na-het'
    else:
      mmpair = 'na-water'

  counts[pairtype] += 1
  mmpair_counts[pairtype][mmpair] += 1
  dotcounts[pairtype] += int(dotcount)
  mmpair_dotcounts[pairtype][mmpair] += int(dotcount)
  if mingap <= -0.5:
    strongcounts[pairtype] += 1
    mmpair_strongcounts[pairtype][mmpair] += 1
  elif mingap > -0.1:
    weakcounts[pairtype] += 1
    mmpair_weakcounts[pairtype][mmpair] += 1
  else:
    midcounts[pairtype] += 1
    mmpair_midcounts[pairtype][mmpair] += 1

mm_res_count = residue_counts[PROTEIN] + residue_counts[NA]

#print_order = ['mcmc','mcsc','scsc','het','water']
#sys.stdout.write("\nCount, all H bonds\n")
#for pairtype in print_order:
#  pass
  #sys.stdout.write("%s: %i, %.2f per residue\n" % (pairtype, counts[pairtype]))
#sys.stdout.write("\nDotcounts\n")
#for pairtype in print_order:
#  sys.stdout.write("%s: %i\n" % (pairtype, dotcounts[pairtype]))
#sys.stdout.write("\nBonds 0.5 mingap or better\n")
#for pairtype in print_order:
#  sys.stdout.write("%s: %i\n" % (pairtype, strongcounts[pairtype]))
#sys.stdout.write("\nBonds between 0.5 and 0.1 mingap\n")
#for pairtype in print_order:
#  sys.stdout.write("%s: %i\n" % (pairtype, midcounts[pairtype]))
#sys.stdout.write("\nBonds worse than 0.1 mingap\n")
#for pairtype in print_order:
#  sys.stdout.write("%s: %i\n" % (pairtype, weakcounts[pairtype]))

#sys.stdout.write()

sys.stdout.write("\n---------------META---------------\n")
#-----------------------META---------------------------------------------
sys.stdout.write("\nFilename: "+os.path.basename(pdbfile)+"\n")
resolution = find_resolution(pdbfile)
if resolution is not None:
  sys.stdout.write("Resolution: "+resolution+"\n")
else:
  sys.stdout.write("No X-ray Resolution found\n")

sys.stdout.write("\nResidue counts\n")
sys.stdout.write("Macromolecule: %i\n" % (residue_counts[PROTEIN]+residue_counts[NA]))
sys.stdout.write("Protein: %i\n" % residue_counts[PROTEIN])
sys.stdout.write("Nucleic acid: %i\n" % residue_counts[NA])
sys.stdout.write("Hets: %i\n" % residue_counts[HET])
sys.stdout.write("Waters: %i\n" % residue_counts[WATER])
#-------------------------------------------------------------------------------

sys.stdout.write("\n---------------BOND COUNTS---------------\n")

#-----------------------BOND COUNTS-----------------------------
sys.stdout.write("\nMacromolecule H bonds\n")
sys.stdout.write("mcmc: %i, %.2f per residue\n" % (counts['mcmc'], (counts['mcmc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("mcsc: %i, %.2f per residue\n" % (counts['mcsc'], (counts['mcsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("scsc: %i, %.2f per residue\n" % (counts['scsc'], (counts['scsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all mc: %i, %.2f per residue\n" % ((counts['mcmc']+counts['mcsc']), ((counts['mcmc']+counts['mcsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all: %i, %.2f per residue\n" % ((counts['mcmc']+counts['mcsc']+counts['scsc']), ((counts['mcmc']+counts['mcsc']+counts['scsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[PROTEIN]:
  sys.stdout.write("\nProtein H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_counts['mcmc']['protein'], (mmpair_counts['mcmc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_counts['mcsc']['protein'], (mmpair_counts['mcsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_counts['scsc']['protein'], (mmpair_counts['scsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_counts['mcmc']['protein']+mmpair_counts['mcsc']['protein']), ((mmpair_counts['mcmc']['protein']+mmpair_counts['mcsc']['protein'])/residue_counts[PROTEIN])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_counts['mcmc']['protein']+mmpair_counts['mcsc']['protein']+mmpair_counts['scsc']['protein']), ((mmpair_counts['mcmc']['protein']+mmpair_counts['mcsc']['protein']+mmpair_counts['scsc']['protein'])/residue_counts[PROTEIN])))

if residue_counts[NA]:
  sys.stdout.write("\nNucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_counts['mcmc']['na'], (mmpair_counts['mcmc']['na']/residue_counts[NA])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_counts['mcsc']['na'], (mmpair_counts['mcsc']['na']/residue_counts[NA])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_counts['scsc']['na'], (mmpair_counts['scsc']['na']/residue_counts[NA])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_counts['mcmc']['na']+mmpair_counts['mcsc']['na']), ((mmpair_counts['mcmc']['na']+mmpair_counts['mcsc']['na'])/residue_counts[NA])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_counts['mcmc']['na']+mmpair_counts['mcsc']['na']+mmpair_counts['scsc']['na']), ((mmpair_counts['mcmc']['na']+mmpair_counts['mcsc']['na']+mmpair_counts['scsc']['na'])/residue_counts[NA])))

if residue_counts[PROTEIN] and residue_counts[NA]:
  sys.stdout.write("\nProtein-Nucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_counts['mcmc']['protein-na'], (mmpair_counts['mcmc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_counts['mcsc']['protein-na'], (mmpair_counts['mcsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_counts['scsc']['protein-na'], (mmpair_counts['scsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_counts['mcmc']['protein-na']+mmpair_counts['mcsc']['protein-na']), ((mmpair_counts['mcmc']['protein-na']+mmpair_counts['mcsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_counts['mcmc']['protein-na']+mmpair_counts['mcsc']['protein-na']+mmpair_counts['scsc']['protein-na']), ((mmpair_counts['mcmc']['protein-na']+mmpair_counts['mcsc']['protein-na']+mmpair_counts['scsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[HET]:
  sys.stdout.write("\nHet H bonds\n")
  sys.stdout.write("het total: %i, %.2f per het\n" % (counts['het'], (counts['het']/(residue_counts[HET]))))
  sys.stdout.write("het-protein: %i, %.2f per het\n" % (mmpair_counts['het']['protein-het'], (mmpair_counts['het']['protein-het']/residue_counts[HET])))
  sys.stdout.write("het-na: %i, %.2f per het\n" % (mmpair_counts['het']['na-het'], (mmpair_counts['het']['na-het']/residue_counts[HET])))

if residue_counts[WATER]:
  sys.stdout.write("\nWater H bonds\n")
  sys.stdout.write("water total: %i, %.2f per HOH\n" % (counts['water'], (counts['water']/(residue_counts[WATER]))))
  sys.stdout.write("water-protein: %i, %.2f per HOH\n" % (mmpair_counts['water']['protein-water'], (mmpair_counts['water']['protein-water']/residue_counts[WATER])))
  sys.stdout.write("water-na: %i, %.2f per HOH\n" % (mmpair_counts['water']['na-water'], (mmpair_counts['water']['na-water']/residue_counts[WATER])))
#-------------------------------------------------------------------------------

sys.stdout.write("\n---------------DOT COUNTS---------------\n")

#---------------------------DOT COUNTS----------------------------------------
sys.stdout.write("\nMacromolecule Dotcounts\n")
sys.stdout.write("mcmc: %i, %.2f per residue\n" % (dotcounts['mcmc'], (dotcounts['mcmc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("mcsc: %i, %.2f per residue\n" % (dotcounts['mcsc'], (dotcounts['mcsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("scsc: %i, %.2f per residue\n" % (dotcounts['scsc'], (dotcounts['scsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all mc: %i, %.2f per residue\n" % ((dotcounts['mcmc']+dotcounts['mcsc']), ((dotcounts['mcmc']+dotcounts['mcsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all: %i, %.2f per residue\n" % ((dotcounts['mcmc']+dotcounts['mcsc']+dotcounts['scsc']), ((dotcounts['mcmc']+dotcounts['mcsc']+dotcounts['scsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[PROTEIN]:
  sys.stdout.write("\nProtein Dotcounts\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_dotcounts['mcmc']['protein'], (mmpair_dotcounts['mcmc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_dotcounts['mcsc']['protein'], (mmpair_dotcounts['mcsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_dotcounts['scsc']['protein'], (mmpair_dotcounts['scsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_dotcounts['mcmc']['protein']+mmpair_dotcounts['mcsc']['protein']), ((mmpair_dotcounts['mcmc']['protein']+mmpair_dotcounts['mcsc']['protein'])/residue_counts[PROTEIN])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_dotcounts['mcmc']['protein']+mmpair_dotcounts['mcsc']['protein']+mmpair_dotcounts['scsc']['protein']), ((mmpair_dotcounts['mcmc']['protein']+mmpair_dotcounts['mcsc']['protein']+mmpair_dotcounts['scsc']['protein'])/residue_counts[PROTEIN])))

if residue_counts[NA]:
  sys.stdout.write("\nNucleic Acid Dotcounts\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_dotcounts['mcmc']['na'], (mmpair_dotcounts['mcmc']['na']/residue_counts[NA])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_dotcounts['mcsc']['na'], (mmpair_dotcounts['mcsc']['na']/residue_counts[NA])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_dotcounts['scsc']['na'], (mmpair_dotcounts['scsc']['na']/residue_counts[NA])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_dotcounts['mcmc']['na']+mmpair_dotcounts['mcsc']['na']), ((mmpair_dotcounts['mcmc']['na']+mmpair_dotcounts['mcsc']['na'])/residue_counts[NA])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_dotcounts['mcmc']['na']+mmpair_dotcounts['mcsc']['na']+mmpair_dotcounts['scsc']['na']), ((mmpair_dotcounts['mcmc']['na']+mmpair_dotcounts['mcsc']['na']+mmpair_dotcounts['scsc']['na'])/residue_counts[NA])))

if residue_counts[PROTEIN] and residue_counts[NA]:
  sys.stdout.write("\nProtein-Nucleic Acid Dotcounts\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_dotcounts['mcmc']['protein-na'], (mmpair_dotcounts['mcmc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_dotcounts['mcsc']['protein-na'], (mmpair_dotcounts['mcsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_dotcounts['scsc']['protein-na'], (mmpair_dotcounts['scsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_dotcounts['mcmc']['protein-na']+mmpair_dotcounts['mcsc']['protein-na']), ((mmpair_dotcounts['mcmc']['protein-na']+mmpair_dotcounts['mcsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_dotcounts['mcmc']['protein-na']+mmpair_dotcounts['mcsc']['protein-na']+mmpair_dotcounts['scsc']['protein-na']), ((mmpair_dotcounts['mcmc']['protein-na']+mmpair_dotcounts['mcsc']['protein-na']+mmpair_dotcounts['scsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[HET]:
  sys.stdout.write("\nHet Dotcounts\n")
  sys.stdout.write("het total: %i, %.2f per het\n" % (dotcounts['het'], (dotcounts['het']/(residue_counts[HET]))))
  sys.stdout.write("het-protein: %i, %.2f per het\n" % (mmpair_dotcounts['het']['protein-het'], (mmpair_dotcounts['het']['protein-het']/residue_counts[HET])))
  sys.stdout.write("het-na: %i, %.2f per het\n" % (mmpair_dotcounts['het']['na-het'], (mmpair_dotcounts['het']['na-het']/residue_counts[HET])))

if residue_counts[WATER]:
  sys.stdout.write("\nWater Dotcounts\n")
  sys.stdout.write("water total: %i, %.2f per HOH\n" % (dotcounts['water'], (dotcounts['water']/(residue_counts[WATER]))))
  sys.stdout.write("water-protein: %i, %.2f per HOH\n" % (mmpair_dotcounts['water']['protein-water'], (mmpair_dotcounts['water']['protein-water']/residue_counts[WATER])))
  sys.stdout.write("water-na: %i, %.2f per HOH\n" % (mmpair_dotcounts['water']['na-water'], (mmpair_dotcounts['water']['na-water']/residue_counts[WATER])))
#-------------------------------------------------------------------------------

sys.stdout.write("\n----------STRONG BOND COUNTS (>=0.5A OVERLAP)-----------\n")

#-----------------------STRONG BOND COUNTS-----------------------------
sys.stdout.write("\nMacromolecule H bonds\n")
sys.stdout.write("mcmc: %i, %.2f per residue\n" % (strongcounts['mcmc'], (strongcounts['mcmc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("mcsc: %i, %.2f per residue\n" % (strongcounts['mcsc'], (strongcounts['mcsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("scsc: %i, %.2f per residue\n" % (strongcounts['scsc'], (strongcounts['scsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all mc: %i, %.2f per residue\n" % ((strongcounts['mcmc']+strongcounts['mcsc']), ((strongcounts['mcmc']+strongcounts['mcsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all: %i, %.2f per residue\n" % ((strongcounts['mcmc']+strongcounts['mcsc']+strongcounts['scsc']), ((strongcounts['mcmc']+strongcounts['mcsc']+strongcounts['scsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[PROTEIN]:
  sys.stdout.write("\nProtein H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_strongcounts['mcmc']['protein'], (mmpair_strongcounts['mcmc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_strongcounts['mcsc']['protein'], (mmpair_strongcounts['mcsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_strongcounts['scsc']['protein'], (mmpair_strongcounts['scsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_strongcounts['mcmc']['protein']+mmpair_strongcounts['mcsc']['protein']), ((mmpair_strongcounts['mcmc']['protein']+mmpair_strongcounts['mcsc']['protein'])/residue_counts[PROTEIN])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_strongcounts['mcmc']['protein']+mmpair_strongcounts['mcsc']['protein']+mmpair_strongcounts['scsc']['protein']), ((mmpair_strongcounts['mcmc']['protein']+mmpair_strongcounts['mcsc']['protein']+mmpair_strongcounts['scsc']['protein'])/residue_counts[PROTEIN])))

if residue_counts[NA]:
  sys.stdout.write("\nNucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_strongcounts['mcmc']['na'], (mmpair_strongcounts['mcmc']['na']/residue_counts[NA])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_strongcounts['mcsc']['na'], (mmpair_strongcounts['mcsc']['na']/residue_counts[NA])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_strongcounts['scsc']['na'], (mmpair_strongcounts['scsc']['na']/residue_counts[NA])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_strongcounts['mcmc']['na']+mmpair_strongcounts['mcsc']['na']), ((mmpair_strongcounts['mcmc']['na']+mmpair_strongcounts['mcsc']['na'])/residue_counts[NA])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_strongcounts['mcmc']['na']+mmpair_strongcounts['mcsc']['na']+mmpair_strongcounts['scsc']['na']), ((mmpair_strongcounts['mcmc']['na']+mmpair_strongcounts['mcsc']['na']+mmpair_strongcounts['scsc']['na'])/residue_counts[NA])))

if residue_counts[PROTEIN] and residue_counts[NA]:
  sys.stdout.write("\nProtein-Nucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_strongcounts['mcmc']['protein-na'], (mmpair_strongcounts['mcmc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_strongcounts['mcsc']['protein-na'], (mmpair_strongcounts['mcsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_strongcounts['scsc']['protein-na'], (mmpair_strongcounts['scsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_strongcounts['mcmc']['protein-na']+mmpair_strongcounts['mcsc']['protein-na']), ((mmpair_strongcounts['mcmc']['protein-na']+mmpair_strongcounts['mcsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_strongcounts['mcmc']['protein-na']+mmpair_strongcounts['mcsc']['protein-na']+mmpair_strongcounts['scsc']['protein-na']), ((mmpair_strongcounts['mcmc']['protein-na']+mmpair_strongcounts['mcsc']['protein-na']+mmpair_strongcounts['scsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[HET]:
  sys.stdout.write("\nHet H bonds\n")
  sys.stdout.write("het total: %i, %.2f per het\n" % (strongcounts['het'], (strongcounts['het']/(residue_counts[HET]))))
  sys.stdout.write("het-protein: %i, %.2f per het\n" % (mmpair_strongcounts['het']['protein-het'], (mmpair_strongcounts['het']['protein-het']/residue_counts[HET])))
  sys.stdout.write("het-na: %i, %.2f per het\n" % (mmpair_strongcounts['het']['na-het'], (mmpair_strongcounts['het']['na-het']/residue_counts[HET])))

if residue_counts[WATER]:
  sys.stdout.write("\nWater H bonds\n")
  sys.stdout.write("water total: %i, %.2f per HOH\n" % (strongcounts['water'], (strongcounts['water']/(residue_counts[WATER]))))
  sys.stdout.write("water-protein: %i, %.2f per HOH\n" % (mmpair_strongcounts['water']['protein-water'], (mmpair_strongcounts['water']['protein-water']/residue_counts[WATER])))
  sys.stdout.write("water-na: %i, %.2f per HOH\n" % (mmpair_strongcounts['water']['na-water'], (mmpair_strongcounts['water']['na-water']/residue_counts[WATER])))
#-------------------------------------------------------------------------------

sys.stdout.write("\n----------MID BOND COUNTS----------\n")

#-----------------------MID BOND COUNTS-----------------------------
sys.stdout.write("\nMacromolecule H bonds\n")
sys.stdout.write("mcmc: %i, %.2f per residue\n" % (midcounts['mcmc'], (midcounts['mcmc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("mcsc: %i, %.2f per residue\n" % (midcounts['mcsc'], (midcounts['mcsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("scsc: %i, %.2f per residue\n" % (midcounts['scsc'], (midcounts['scsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all mc: %i, %.2f per residue\n" % ((midcounts['mcmc']+midcounts['mcsc']), ((midcounts['mcmc']+midcounts['mcsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all: %i, %.2f per residue\n" % ((midcounts['mcmc']+midcounts['mcsc']+midcounts['scsc']), ((midcounts['mcmc']+midcounts['mcsc']+midcounts['scsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[PROTEIN]:
  sys.stdout.write("\nProtein H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_midcounts['mcmc']['protein'], (mmpair_midcounts['mcmc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_midcounts['mcsc']['protein'], (mmpair_midcounts['mcsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_midcounts['scsc']['protein'], (mmpair_midcounts['scsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_midcounts['mcmc']['protein']+mmpair_midcounts['mcsc']['protein']), ((mmpair_midcounts['mcmc']['protein']+mmpair_midcounts['mcsc']['protein'])/residue_counts[PROTEIN])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_midcounts['mcmc']['protein']+mmpair_midcounts['mcsc']['protein']+mmpair_midcounts['scsc']['protein']), ((mmpair_midcounts['mcmc']['protein']+mmpair_midcounts['mcsc']['protein']+mmpair_midcounts['scsc']['protein'])/residue_counts[PROTEIN])))

if residue_counts[NA]:
  sys.stdout.write("\nNucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_midcounts['mcmc']['na'], (mmpair_midcounts['mcmc']['na']/residue_counts[NA])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_midcounts['mcsc']['na'], (mmpair_midcounts['mcsc']['na']/residue_counts[NA])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_midcounts['scsc']['na'], (mmpair_midcounts['scsc']['na']/residue_counts[NA])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_midcounts['mcmc']['na']+mmpair_midcounts['mcsc']['na']), ((mmpair_midcounts['mcmc']['na']+mmpair_midcounts['mcsc']['na'])/residue_counts[NA])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_midcounts['mcmc']['na']+mmpair_midcounts['mcsc']['na']+mmpair_midcounts['scsc']['na']), ((mmpair_midcounts['mcmc']['na']+mmpair_midcounts['mcsc']['na']+mmpair_midcounts['scsc']['na'])/residue_counts[NA])))

if residue_counts[PROTEIN] and residue_counts[NA]:
  sys.stdout.write("\nProtein-Nucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_midcounts['mcmc']['protein-na'], (mmpair_midcounts['mcmc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_midcounts['mcsc']['protein-na'], (mmpair_midcounts['mcsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_midcounts['scsc']['protein-na'], (mmpair_midcounts['scsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_midcounts['mcmc']['protein-na']+mmpair_midcounts['mcsc']['protein-na']), ((mmpair_midcounts['mcmc']['protein-na']+mmpair_midcounts['mcsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_midcounts['mcmc']['protein-na']+mmpair_midcounts['mcsc']['protein-na']+mmpair_midcounts['scsc']['protein-na']), ((mmpair_midcounts['mcmc']['protein-na']+mmpair_midcounts['mcsc']['protein-na']+mmpair_midcounts['scsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[HET]:
  sys.stdout.write("\nHet H bonds\n")
  sys.stdout.write("het total: %i, %.2f per het\n" % (midcounts['het'], (midcounts['het']/(residue_counts[HET]))))
  sys.stdout.write("het-protein: %i, %.2f per het\n" % (mmpair_midcounts['het']['protein-het'], (mmpair_midcounts['het']['protein-het']/residue_counts[HET])))
  sys.stdout.write("het-na: %i, %.2f per het\n" % (mmpair_midcounts['het']['na-het'], (mmpair_midcounts['het']['na-het']/residue_counts[HET])))

if residue_counts[WATER]:
  sys.stdout.write("\nWater H bonds\n")
  sys.stdout.write("water total: %i, %.2f per HOH\n" % (midcounts['water'], (midcounts['water']/(residue_counts[WATER]))))
  sys.stdout.write("water-protein: %i, %.2f per HOH\n" % (mmpair_midcounts['water']['protein-water'], (mmpair_midcounts['water']['protein-water']/residue_counts[WATER])))
  sys.stdout.write("water-na: %i, %.2f per HOH\n" % (mmpair_midcounts['water']['na-water'], (mmpair_midcounts['water']['na-water']/residue_counts[WATER])))
#-------------------------------------------------------------------------------

sys.stdout.write("\n----------WEAK BOND COUNTS (<0.1A OVERLAP)----------\n")

#-----------------------WEAK BOND COUNTS-----------------------------
sys.stdout.write("\nMacromolecule H bonds\n")
sys.stdout.write("mcmc: %i, %.2f per residue\n" % (weakcounts['mcmc'], (weakcounts['mcmc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("mcsc: %i, %.2f per residue\n" % (weakcounts['mcsc'], (weakcounts['mcsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("scsc: %i, %.2f per residue\n" % (weakcounts['scsc'], (weakcounts['scsc']/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all mc: %i, %.2f per residue\n" % ((weakcounts['mcmc']+weakcounts['mcsc']), ((weakcounts['mcmc']+weakcounts['mcsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
sys.stdout.write("all: %i, %.2f per residue\n" % ((weakcounts['mcmc']+weakcounts['mcsc']+weakcounts['scsc']), ((weakcounts['mcmc']+weakcounts['mcsc']+weakcounts['scsc'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[PROTEIN]:
  sys.stdout.write("\nProtein H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_weakcounts['mcmc']['protein'], (mmpair_weakcounts['mcmc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_weakcounts['mcsc']['protein'], (mmpair_weakcounts['mcsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_weakcounts['scsc']['protein'], (mmpair_weakcounts['scsc']['protein']/residue_counts[PROTEIN])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_weakcounts['mcmc']['protein']+mmpair_weakcounts['mcsc']['protein']), ((mmpair_weakcounts['mcmc']['protein']+mmpair_weakcounts['mcsc']['protein'])/residue_counts[PROTEIN])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_weakcounts['mcmc']['protein']+mmpair_weakcounts['mcsc']['protein']+mmpair_weakcounts['scsc']['protein']), ((mmpair_weakcounts['mcmc']['protein']+mmpair_weakcounts['mcsc']['protein']+mmpair_weakcounts['scsc']['protein'])/residue_counts[PROTEIN])))

if residue_counts[NA]:
  sys.stdout.write("\nNucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_weakcounts['mcmc']['na'], (mmpair_weakcounts['mcmc']['na']/residue_counts[NA])))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_weakcounts['mcsc']['na'], (mmpair_weakcounts['mcsc']['na']/residue_counts[NA])))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_weakcounts['scsc']['na'], (mmpair_weakcounts['scsc']['na']/residue_counts[NA])))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_weakcounts['mcmc']['na']+mmpair_weakcounts['mcsc']['na']), ((mmpair_weakcounts['mcmc']['na']+mmpair_weakcounts['mcsc']['na'])/residue_counts[NA])))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_weakcounts['mcmc']['na']+mmpair_weakcounts['mcsc']['na']+mmpair_weakcounts['scsc']['na']), ((mmpair_weakcounts['mcmc']['na']+mmpair_weakcounts['mcsc']['na']+mmpair_weakcounts['scsc']['na'])/residue_counts[NA])))

if residue_counts[PROTEIN] and residue_counts[NA]:
  sys.stdout.write("\nProtein-Nucleic Acid H bonds\n")
  sys.stdout.write("mcmc: %i, %.2f per residue\n" % (mmpair_weakcounts['mcmc']['protein-na'], (mmpair_weakcounts['mcmc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("mcsc: %i, %.2f per residue\n" % (mmpair_weakcounts['mcsc']['protein-na'], (mmpair_weakcounts['mcsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("scsc: %i, %.2f per residue\n" % (mmpair_weakcounts['scsc']['protein-na'], (mmpair_weakcounts['scsc']['protein-na']/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all mc: %i, %.2f per residue\n" % ((mmpair_weakcounts['mcmc']['protein-na']+mmpair_weakcounts['mcsc']['protein-na']), ((mmpair_weakcounts['mcmc']['protein-na']+mmpair_weakcounts['mcsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))
  sys.stdout.write("all: %i, %.2f per residue\n" % ((mmpair_weakcounts['mcmc']['protein-na']+mmpair_weakcounts['mcsc']['protein-na']+mmpair_weakcounts['scsc']['protein-na']), ((mmpair_weakcounts['mcmc']['protein-na']+mmpair_weakcounts['mcsc']['protein-na']+mmpair_weakcounts['scsc']['protein-na'])/(residue_counts[PROTEIN]+residue_counts[NA]))))

if residue_counts[HET]:
  sys.stdout.write("\nHet H bonds\n")
  sys.stdout.write("het total: %i, %.2f per het\n" % (weakcounts['het'], (weakcounts['het']/(residue_counts[HET]))))
  sys.stdout.write("het-protein: %i, %.2f per het\n" % (mmpair_weakcounts['het']['protein-het'], (mmpair_weakcounts['het']['protein-het']/residue_counts[HET])))
  sys.stdout.write("het-na: %i, %.2f per het\n" % (mmpair_weakcounts['het']['na-het'], (mmpair_weakcounts['het']['na-het']/residue_counts[HET])))

if residue_counts[WATER]:
  sys.stdout.write("\nWater H bonds\n")
  sys.stdout.write("water total: %i, %.2f per HOH\n" % (weakcounts['water'], (weakcounts['water']/(residue_counts[WATER]))))
  sys.stdout.write("water-protein: %i, %.2f per HOH\n" % (mmpair_weakcounts['water']['protein-water'], (mmpair_weakcounts['water']['protein-water']/residue_counts[WATER])))
  sys.stdout.write("water-na: %i, %.2f per HOH\n" % (mmpair_weakcounts['water']['na-water'], (mmpair_weakcounts['water']['na-water']/residue_counts[WATER])))
#-------------------------------------------------------------------------------

sys.stdout.write("\n")

#ToDo
#count donor/acceptors, N and O, or fixed per sidechain?  how to alts?
#helix and sheet
