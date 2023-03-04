import sys

def pdic2p(pdic):
  """
  Converts the dictionary {'XX' : 1.2 , 'XY' : -1.6 ...} to an ordered list}
  """
  keys = list(pdic.keys())
  keys.sort()
  sorted_pdic = {i: pdic[i] for i in keys}

  return list(sorted_pdic.values())

if len(sys.argv) != 2:
  print("Usage: python3 {} <qhem_output>".format(sys.argv[0]))
  sys.exit()

# Extract information
mp = []
xyz = []
ats = []
with open(sys.argv[1],'r') as ff:
  while True:
    line = ff.readline()
    if line == "":
      break

    if "  Standard Nuclear" in line:
      ff.readline()
      ff.readline()
      line = ff.readline()
      while not "----------" in line:
        ls = line.strip().split()
        xyz.append([float(i) for i in ls[2:]])
        ats.append(ls[1])
        line = ff.readline()

    if "  The traceless molecular multipole moments" in line:
      ff.readline()  
      ff.readline()
      ff.readline()
      ff.readline()

      # Assuming charge and dipoles are always there
      chg = float(ff.readline().strip().split()[0])
      mp.append([chg])
      ff.readline()

      # Go for dipoles
      dip_line = ff.readline().strip().split()
      p = [float(dip_line[1]),float(dip_line[3])]
      dip_line = ff.readline().strip().split()
      p.append(float(dip_line[1]))

      pdic = {}
      flag = False
      
      line = ff.readline()
      # Go for the rest of multipoles
      while not "-------" in line:
        if "Debye-" in line:
          if flag:
            p = list(pdic.values())
          flag = True
          mp.append(p)
          p = []
          pdic = {}
          line = ff.readline()

        sl = line.strip().split()
        for i in range(0,len(sl),2):
          pdic[sl[i]] = float(sl[i+1])
        line = ff.readline()

      if flag:
        p = list(pdic.values())
      mp.append(p)


mp_converted_eA = [[chg]]
for i in range(1,len(mp)):
  mp_converted_eA.append([mp[i][j] * 0.2081943 for j in range(len(mp[i]))])

with open('training_set_frame.xyz','w') as ff:
  ff.write("{}\n".format(str(len(ats))))
  ff.write("{}\n".format(str(mp_converted_eA)))
  for i in range(len(ats)):
    ff.write("{0:5s} {1:20.8f} {2:20.8f} {3:20.8f}\n".format(ats[i],xyz[i][0],xyz[i][1],xyz[i][2]))

