import os
import sys
import numpy as np

def get_indices(prompt,max_index):
    while True:
        user_input=input(prompt).strip()
        if not user_input:
            print("Error: no indices provided. Please enter at least one integer.")
            continue
        try:
            indices=list(map(int,user_input.split()))
        except ValueError:
            print("Error: please enter only integer numbers separated by spaces.")
            continue
        if len(indices)!=len(set(indices)):
            print("Error: duplicate indices detected. Please enter unique indices.")
            continue
        if any(i<0 or i>=max_index for i in indices):
            print(f"Error: indices must be between 0 and {max_index-1}.")
            continue
        return indices

### Definition of the atomic mass dictionary

element_values={
    "C": 12.01000,
    "H": 1.00800,
    "O": 16.00000,
    "N": 14.01000,
    "S": 32.06000,
    "P": 30.97000,
    "Na": 22.98980 
}

### Warning about linear or triatomic molecules
print("\n************************************************************")
print("WARNING")
print("The PMM-QMMM method may produce unreliable results for:")
print(" • Triatomic molecules")
print(" • Perfectly linear molecules")
print("Please verify that your system is suitable before proceeding.")
print("************************************************************\n")

### Input and output files

pdb=input("Enter the PDB file name:\n").strip()
if not os.path.exists(pdb):
    print(f"Error: Missing file {pdb}")
    sys.exit(1)
prms=input("Enter the ORCA forcefield file (*prms):\n").strip()
if not os.path.exists(prms):
    print(f"Error: Missing file {prms}")
    sys.exit(1)
ch=input("Enter the charge of the Quantum Center (QC):\n").strip()
try:
    ch=int(ch)
except ValueError:
    print("The entered value is not an integer!")
    sys.exit(1)
sc="MD_charges.txt" # MD simulation charges file
qcm="QC_indexes.txt" # Quantum Center (QC) index file with masses
ef="env_indexes.txt" # Environment index file
fqc="file_QC.txt" # QM geometry file
pc=input("Enter the ORCA point-charges file:\n").strip()
if not os.path.exists(pc):
    print(f"Error: Missing file {pc}")
    sys.exit(1)
xtc=input("Enter the MD trajectory file in xtc format:\n")
if not os.path.exists(xtc):
    print(f"Error: Missing file {xtc}")
    sys.exit(1)
ff=input("Enter the initial frame (0 for the first frame):\n")
try:
    ff=int(ff)
except ValueError:
    print("The entered value is not an integer!")
    sys.exit(1)
ff=str(ff)
lf=input("Enter the last frame:\n")
try:
    lf=int(lf)
except ValueError:
    print("The entered value is not an integer!")
    sys.exit(1)
lf=str(lf)

with open(prms,'r') as f:
    lps=f.readlines()
na=int(lps[3].split()[0]) # Number of atoms
lps=lps[4:]
ind1=get_indices("Enter the indices of the atoms in the Quantum Center (QC) (excluding the saturating atom), starting from 0:\n",na)
nqc=len(ind1) # Number of atoms in the Quantum Center

### Generation of the parameter files

cha=[0 for _ in range(na)]
chb=[0 for _ in range(na)]
elem=[0 for _ in range(na)]
mass=[0 for _ in range(na)]
for k in range(na):
    cha[k]=float(lps[k].split()[2])
    chb[k]=float(lps[k].split()[2])
    elem[k]=lps[k].split()[1]
    try:
        mass[k]=float(element_values[elem[k]])
    except KeyError:
        print(f"Warning: Element '{elem[k]}' not found in 'element_values'.")
        print("Please add its atomic mass in the 'element_values' dictionary at the beginning of the program.")
        sys.exit(1)

q=input("Are there any covalent bonds between the Quantum Center (QC) and atoms outside the QC? (y/n)\n")
if q=="y":
    exc=0 # Excess charge in the QC from topology parameters
    for t in ind1:
        exc=exc+cha[t]
        chb[t]=0.0
    ind2=get_indices("Enter the indices of the atoms directly bonded to the Quantum Center (QC), starting from 0:\n",na)
    for t2 in ind2:
        exc=exc+cha[t]
        chb[t2]=0.0
    exc=exc-ch
    ind3=get_indices("Enter the indices of atoms belonging to the Quantum Center (QC) fragment but outside the QC itself, starting from 0:\n",na)
    ind=[x for x in ind3 if x not in ind2]
    for t in ind:
        chb[t]=cha[t]-(exc/len(ind))
elif q=="n":
    for t in ind1:
        chb[t]=0.0
else:
    print("Invalid selection. Please enter y or n.")
out1=open(sc,'w')
out2=open(qcm,'w')
out3=open(ef,'w')
out2.write(str(len(ind1))+"\n")
out3.write(str(na-len(ind1))+"\n")
for k in range(na):
    out1.write("{:>7d}{:>12.6f}".format(k+1,chb[k])+"\n")
    if k in ind1:
        out2.write("{:>7d}{:>10.5f}".format(k+1,mass[k])+"\n")
    else:
        out3.write("{:>7d}".format(k+1)+"\n")
out1.close()
out2.close()
out3.close()
ll=open(pdb,'r').readlines()
out=open(fqc,'w')
out.write(str(len(ind1))+"\n")
for l in ll:
    if l.startswith("ATOM") or l.startswith("HETATM"):
        ai=int(l[6:11].strip())-1
        if ai in ind1:
            x=float(l[30:38].strip())
            y=float(l[38:46].strip())
            z=float(l[46:54].strip())
            e=elem[ai]
            m=float(element_values[e])#
            out.write("{:>10.5f}{:>15.3f}{:>15.3f}{:>15.3f}".format(m,x,y,z)+"\n")
q=input("Do you want to perform the fit on a subregion of the Quantum Center (QC)? (y/n)\n")
if q=="y":
    ind5=get_indices("Enter the indices of the atoms in the selected subregion of the Quantum Center (QC), starting from 0:\n",na)
elif q=="n":
    ind5=ind1
else:
    print("Invalid selection. Please enter y or n.")
out.write(str(len(ind5))+"\n")
for l in ll:
    if l.startswith("ATOM") or l.startswith("HETATM"):
        ai=int(l[6:11].strip())-1
        if ai in ind5:
            x=float(l[30:38].strip())
            y=float(l[38:46].strip())
            z=float(l[46:54].strip())
            e=elem[ai]
            m=float(element_values[e])
            out.write("{:>10.5f}{:>15.3f}{:>15.3f}{:>15.3f}".format(m,x,y,z)+"\n")
out.close()

### Generation of the input file for PMM-QMMM

n=input("Enter the number of electronic excited states (enter 0 for ground state only):\n")
try:
    n=int(n)
except ValueError:
    print("The entered value is not an integer!")
    sys.exit(1)
n=n+1 # The total number of electronic states is n + 1, including the ground state
mux=[[0 for _ in range(n)] for _ in range(n)] # Matrix of the x-component of the electric dipole moment
muy=[[0 for _ in range(n)] for _ in range(n)] # Matrix of the y-component of the electric dipole moment
muz=[[0 for _ in range(n)] for _ in range(n)] # Matrix of the z-component of the electric dipole moment
ene=[0 for _ in range(n)] # vector of energies
for i in range(n):
    file="S"+str(i)+"/ex_state.out"
    if not os.path.exists(file):
        print(f"Error: Missing file {file}")
        sys.exit(1)
    lines=open(file,'r').readlines()
    j=1
    check1=True
    check2=True
    while j<len(lines) and check2:
        line=lines[-j]
        if "Total Dipole Moment" in line and check1:
            mux[i][i]=float(line.split()[-3])
            muy[i][i]=float(line.split()[-2])
            muz[i][i]=float(line.split()[-1])
            check1=False
        elif "FINAL SINGLE POINT ENERGY   " in lines[-j] and i==0:
            ene[i]=float(line.split()[-1]) # Ground-state energy
            check2=False
        j=j+1
# Construction of the electric dipole moment and energy matrices
file="S"+str(n-1)+"/ex_state.out"
lines=open(file,'r').readlines()
i=0
while i<len(lines):
    if "TRANSIENT ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in lines[i]:
        i=i+5
        while True:
            if len(lines[i].split())==0:
                break
            else:
                elem=lines[i].split()
                h=int(elem[0].split('-')[0])
                e=int(elem[2].split('-')[0])
                if (h<n) and (e<n):
                    mux[h][e]=float(elem[-3])
                    muy[h][e]=float(elem[-2])
                    muz[h][e]=float(elem[-1])
                    mux[e][h]=float(elem[-3])
                    muy[e][h]=float(elem[-2])
                    muz[e][h]=float(elem[-1])
                i=i+1
    elif "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in lines[i]:
        i=i+5
        while True:
            if len(lines[i].split())==0:
                break
            else:
                elem=lines[i].split()
                h=int(elem[0].split('-')[0])
                e=int(elem[2].split('-')[0])
                if (e<n):
                    mux[h][e]=float(elem[-3])
                    muy[h][e]=float(elem[-2])
                    muz[h][e]=float(elem[-1])
                    mux[e][h]=float(elem[-3])
                    muy[e][h]=float(elem[-2])
                    muz[e][h]=float(elem[-1])
                    ene[e]=float(elem[3]) # Excitation energies
                i=i+1
    i=i+1

na=str(na)
ch=str(ch)
q=input("Select the type of calculation:\n [0] Perturbed Electronic State Properties\n [1] Perturbed Properties for Spectral Absorption Bands\n")
if q=="0":
    i=int(q)
elif q=="2":
    i=inr(q)
else:
    print("Invalid selection. Please enter 0 or 1.")

with open("input-pmm-QMMM", 'w') as out:
    out.write(f"""ORCA point-charges file
{pc}
QM geometry file
{fqc}
MD trajectory file in xtc format
{xtc}
Initial frame (0 for the first frame)
{ff}
Last frame
{lf}
Number of atoms
{na}
Simulation charges file
{sc}
Quantum Center (QC) index file with masses
{qcm}
Environment index file
{ef}
Charge of the Quantum Center (QC)
{ch}
Number of electronic states (including the ground state)
{n}
Number of atoms in the Quantum Center (QC)
{nqc}
Matrix of the electric dipole moment
""")

    for a in range(n):
        for b in range(n):
            out.write("{:>4d}{:>4d}{:>15.9f}{:>15.9f}{:>15.9f}".format(a,b,mux[a][b],muy[a][b],muz[a][b])+"\n")
    
    out.write("Ground-state energy (Hartree)\n")
    out.write(str(ene[0])+"\n")
    out.write("Excitation energies (eV)\n")
    for a in range(1,n):
        out.write(str(ene[a])+"\n")
    out.write("Calculation type\n")
    out.write(str(i)+"\n")

print("The input files have been successfully generated!")
