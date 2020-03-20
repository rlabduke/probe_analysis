import sys
from libtbx import easy_run

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

#probe_out = easy_run.fully_buffered(probe_command)
#probe_results = probe_out.stdout_lines

class probe_line():
  """parses a line of probe unformatted output into a python class"""
  def __init__(self, line, skipwaters=True):
    """accepts a probe unformatted output line as its input"""
    if not line.strip(): # averts an IndexError problem with empty lines
      self.val=False
      return
    self.val = True #val is used by __nonzero__ or __bool__ to test whether this object contains anything

    x = line.split(':')
    #column counts are different between outputs, so set a new starting index where the columns change
    if len(x) == 19:
      #regular unformatted output, one line per dot
      i = 5
    elif len(x) == 20:
      #condensed output with dotcount
      self.condensed = True
      i = 6
    else:
      sys.stderr.write("probe line with unexpected column count:\n")
      sys.stderr.write(line+"\n")
      sys.stderr.write("%i columns found, expected 19 or 20" % (len(x)))
      sys.exit()

    self.name = x[0]
    self.pattern = x[1]
    self.interactiontype = x[2]
    #if not interactiontype == 'hb': continue  # skip non-h-bonds

    src = x[3]
    self.srcChain = src[0:2].strip()
    self.srcNum = int(src[2:6].strip())
    self.srcIns = src[6:7]  # .strip()
    self.srcResname = src[7:10].strip()
    if self.srcResname == 'HOH' and skipwaters:
      self.val=False
      return
    self.srcAtom = src[11:15]  # .strip()
    self.srcAlt = src[15:16].strip()
    self.srcResid = src[0:10] + src[15:16]

    trg = x[4]
    # going to count dots per bond as a measure of strength instead
    self.trgChain = trg[0:2].strip()
    self.trgNum = int(trg[2:6].strip())
    self.trgNumStr = trg[2:6]
    self.trgIns = trg[6:7]  # .strip()
    self.trgResname = trg[7:10].strip()
    self.trgAtom = trg[11:15]  # .strip()
    self.trgAlt = trg[15:16].strip()
    self.trgResid = trg[0:10] + trg[15:16]

    # name:pat:type:srcAtom:targAtom:*dotcount*:min-gap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval
    if self.condensed:
      self.dotcount = x[5]
    self.mingap = float(x[i])
    self.gap = float(x[i + 1])
    self.spX = float(x[i + 2]) #probably coordinates of the far end of a spike
    self.spX = float(x[i + 3])
    self.spX = float(x[i + 4])
    self.spikeLen = float(x[i + 5])
    self.score = float(x[i + 6])
    self.stype = x[i + 7]
    self.ttype = x[i + 8]
    self.x = float(x[i + 9])
    self.y = float(x[i + 10])
    self.z = float(x[i + 11])
    self.sBval = float(x[i + 12])
    self.tBval = float(x[i + 13])

  #these methods should allow "if probe_line is True:" constructions to test object content
  #val defaults to True, is set to False on certain fail cases.
  def __nonzero__(self): #python2 version
    return self.val
  def __bool__(self): #python3 version
    return self.val

#probecommand = 'phenix.probe -quiet -u -condensed -self -mc -nowaters ALL '+pdbfile
#probe_out = easy_run.fully_buffered(probecommand)
#probefile = probe_out.stdout_lines

def check_probe_error(probe_out):
  if probe_out.return_code != 0:
    raise RuntimeError("Probe crashed - dumping stderr:\n%s" %
                       "\n".join(probe_out.stderr_lines))

def run_probe_dots_molprobity(pathtofile, nuclear=False, hb=True, vdw=True, clash=True):
  #This is how MolProbity runs probe to generate dots for mkinemage markup in the function makeProbeDots()
  ocutval = 10
  options = []
  if hb: options.append("-nohbout")
  if vdw: options.append("-novdwout")
  if clash: options.append("-noclashout")
  if nuclear: options.append("-nuclear")
  options = " ".join(options)

  #-dotmaster adds a "dots" master -- useful when using this kin with Probe remote update
  probecommand = "phenix.probe %s -sepworse -4H -quiet -noticks -nogroup -dotmaster -mc -het -self 'ogt%i' %s" %  (options, ocutval,pathtofile)
  probe_out = easy_run.fully_buffered(probecommand)
  check_probe_error(probe_out)
  probefile = probe_out.stdout_lines
  return probefile
  #MolProbity also contains this version for if a clashlimit is passed. Not sure it's used
  #exec("phenix.probe $options -sepworse -DIVlow$clashLimit -4H -quiet -noticks -nogroup -dotmaster -mc -het -self 'ogt$ocutval' $infile >> $outfile");

def run_probe_clashscore(pathtofile, nuclear=False, condensed=False, ogt=10, atomdump=False):
  #from phenix.clashscore
  probe_command = 'phenix.probe'
  nuclear_flag = ""
  condensed_flag = ""
  if nuclear:
    nuclear_flag = "-nuclear"
  if condensed:
    condensed_flag = "-CON"
  probe_txt = \
    '%s -u -q -mc -het -once -NOVDWOUT %s %s' % (probe_command, condensed_flag, nuclear_flag) +\
      ' "ogt%d not water" "ogt%d" %s' % (ogt, ogt, pathtofile)
  #The -NOVDWOUT probe run above is faster for clashscore to parse,
  # the full_probe_txt version below is for printing to file for coot usage
  full_probe_txt = \
    '%s -u -q -mc -het -once %s' % (probe_command, nuclear_flag) +\
      ' "ogt%d not water" "ogt%d" %s' % (ogt, ogt, pathtofile)
  probe_atom_txt = \
    '%s -q -mc -het -dumpatominfo %s' % (probe_command, nuclear_flag) +\
      ' "ogt%d not water" %s' % (ogt, pathtofile)
  #if blt is not None:
  #  self.probe_atom_b_factor = \
  #    '%s -q -mc -het -dumpatominfo %s' % (probe_command, nuclear_flag) +\
  #      ' "blt%d ogt%d not water" -' % (blt, ogt)
  if atomdump:
    probe_out = easy_run.fully_buffered(probe_atom_txt)
  elif condensed:
    probe_out = easy_run.fully_buffered(probe_txt)
  else:
    probe_out = easy_run.fully_buffered(full_probe_txt)
  check_probe_error(probe_out)
  probefile = probe_out.stdout_lines
  return probefile

def run_probe_undowser(pathtofile, nuclear=False, clashonly=True):
  nuclear_flag = ""
  if nuclear:
    nuclear_flag = "-nuclear"
  if clashonly:
    probe_command = "phenix.probe -u -q -mc -het -con %s -once -wat2wat -stdbonds -onlybadout 'water' 'all' %s" % (nuclear_flag, pathtofile)
  else:
    probe_command = "phenix.probe -u -q -mc -het -con %s -once -wat2wat -stdbonds 'water' 'all' %s" % (nuclear_flag, pathtofile)
  probe_out = easy_run.fully_buffered(probe_command)
  check_probe_error(probe_out)
  probefile = probe_out.stdout_lines
  return probefile
