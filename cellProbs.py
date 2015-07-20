
#
# Calculates the probabilty of finding an ion absorption
# given the existance of another ion absorption


import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt



# Calculates the probability of finding b given a
def probAgivenB(a, b):

    totalCells = float(len(a))
    probA = float( a.count(1) ) / totalCells
    
    countBoth = 0.
    for i in range(0,len(a)):
        
        if a[i]==1 and b[i]==1:
            countBoth += 1.

    probAandB = countBoth / totalCells
    probAgB = probAandB / probA
    
    return probAgB

hi   = []
mgii = []
civ  = []
ovi  = []

f = open('dat.files')

for fileName in f:

    print fileName.strip()
    fN = fileName.strip()
    
    ftmp = open(fN)
    ftmp.readline()

    for line in ftmp:

        l2 = line.split()
        hiTmp   = float(l2[8])
        mgiiTmp = float(l2[9])
        civTmp  = float(l2[10])
        oviTmp  = float(l2[11])
        
        hi.append(hiTmp)
        mgii.append(mgiiTmp)
        civ.append(civTmp)
        ovi.append(oviTmp)


    ftmp.close()

f.close()


# Read in the list of all probbed cells
f = open('probed.files')
probedCells = []
for fileName in f:

    fN = fileName.strip()
    ftmp = open(fN)
    ftmp.readline()

    for line in ftmp:
    
        cellnum = int(line)
        probedCells.append(cellnum)

    ftmp.close()
f.close()


# Calculate the probabiltiy of finding MgII given that HI is found
totalCells = float(len(probedCells))
#totalCells = float(len(hi))
print ''
print 'Total Number of Cells:          ', totalCells
print 'Number of HI Absorbing Cells:   ', hi.count(1)
print 'Number of MgII Absorbing Cells: ', mgii.count(1)
print 'Number of CIV Absorbing Cells:  ', civ.count(1)
print 'Number of OVI Absorbing Cells:  ', ovi.count(1)
print ''

# Find probability of each ion
probHI   = float(hi.count(1))/totalCells
probMgII = float(mgii.count(1))/totalCells
probCIV  = float(civ.count(1))/totalCells
probOVI  = float(ovi.count(1))/totalCells

# Find the number of occurances of HI and MgII
count = 0.
for i in range(0,len(hi)):

    if hi[i]==1 and mgii[i]==1:
        count+=1

probAandB = count/totalCells

probMgIIgivenHI = probAandB / probHI
print probMgIIgivenHI


probHI = float(hi.count(1))/totalCells 
probMgII = float(mgii.count(1))/totalCells
probCIV = float(civ.count(1))/totalCells
probOVI = float(ovi.count(1))/totalCells

fo = open('probs.out', 'w')
s = 'Given \t Percent Probabiltiy of Finding\n'
fo.write(s)
s = '      \t HI \t MgII \t CIV \t OVI \n'
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('HI', probHI*100, probAgivenB(hi,mgii)*100, probAgivenB(hi,civ)*100, probAgivenB(hi,ovi)*100)
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('MgII', probAgivenB(mgii,hi)*100, probMgII*100, probAgivenB(mgii,civ)*100, probAgivenB(mgii,ovi)*100)
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('CIV',  probAgivenB(civ,hi)*100, probAgivenB(civ,mgii)*100, probCIV*100, probAgivenB(civ,ovi)*100)
fo.write(s)
s = '{0:s} \t {1:>.3} \t {2:>.3} \t {3:>.3} \t {4:>.3} \n'.format('OVI', probAgivenB(ovi,hi)*100, probAgivenB(ovi,mgii)*100, probAgivenB(ovi,civ)*100, probOVI*100)
fo.write(s)

fo.write('\n\n\n')


##########
# Calculate the probabilities for finding each ion if HI is NOT found
##########

# Loop through the files and pull out all ions that have absorption in any ion
# but not HI

f = open('dat.files')
mgii_nH, mgii_t, mgii_z, mgii_l = [], [], [], []
civ_nH, civ_t, civ_z, civ_l = [], [], [], []
ovi_nH, ovi_t, ovi_z, ovi_l = [], [], [], []
mgiiCount = 0.
civCount = 0.
oviCount = 0.
for fN in f:
    ftmp = open(fN.strip())
    ftmp.readline()
    for line in ftmp:
        l = line.split()
        hiFlag = int(l[8])
        mgiiFlag = int(l[9])
        civFlag = int(l[10])
        oviFlag = int(l[11])

        
        if mgiiFlag==1 and hiFlag==0:
            mgii_nH.append( pow(10,float(l[4])))
            mgii_t.append(pow(10, float(l[5])))
            mgii_z.append(float(l[6]))
            mgii_l.append(float(l[7]))
            mgiiCount += 1
    
        if civFlag==1 and hiFlag==0:
            civ_nH.append(pow(10, float(l[4])))
            civ_t.append(pow(10, float(l[5])))
            civ_z.append(float(l[6]))
            civ_l.append(float(l[7]))
            civCount += 1

        if oviFlag==1 and hiFlag==0:
            ovi_nH.append(pow(10, float(l[4])))
            ovi_t.append(pow(10, float(l[5])))
            ovi_z.append(float(l[6]))
            ovi_l.append(float(l[7]))
            oviCount += 1

    ftmp.close()
f.close()



# Look at MgII
# Fraction of mgII cells that do not have hi to those that do
mgFrac = mgiiCount / mgii.count(1)
civFrac = civCount / civ.count(1)
oviFrac = oviCount / ovi.count(1)

# Density
mg_nH_mean = np.log10( np.mean( mgii_nH ) )
civ_nH_mean = np.log10( np.mean( civ_nH ) )
ovi_nH_mean = np.log10( np.mean( ovi_nH ) )

mg_nH_min = np.log10( np.min( mgii_nH ) )
civ_nH_min = np.log10( np.min( civ_nH ) )
ovi_nH_min = np.log10( np.min( ovi_nH ) )

mg_nH_max = np.log10( np.max( mgii_nH ) )
civ_nH_max = np.log10( np.max( civ_nH ) )
ovi_nH_max = np.log10( np.max( ovi_nH ) )

mg_nH_std = np.log10( np.std( mgii_nH ) )
civ_nH_std = np.log10( np.std( civ_nH ) )
ovi_nH_std = np.log10( np.std( ovi_nH ) )

mg_nH_med = np.log10( np.median( mgii_nH ) )
civ_nH_med = np.log10( np.median( civ_nH ) )
ovi_nH_med = np.log10( np.median( ovi_nH ) )


# Temperature
mg_t_mean = np.log10( np.mean( mgii_t ) )
civ_t_mean = np.log10( np.mean( civ_t ) )
ovi_t_mean = np.log10( np.mean( ovi_t ) )

mg_t_min = np.log10( np.min( mgii_t ) )
civ_t_min = np.log10( np.min( civ_t ) )
ovi_t_min = np.log10( np.min( ovi_t ) )

mg_t_max = np.log10( np.max( mgii_t ) )
civ_t_max = np.log10( np.max( civ_t ) )
ovi_t_max = np.log10( np.max( ovi_t ) )

mg_t_std = np.log10( np.std( mgii_t ) )
civ_t_std = np.log10( np.std( civ_t ) )
ovi_t_std = np.log10( np.std( ovi_t ) )

mg_t_med = np.log10( np.median( mgii_t ) )
civ_t_med = np.log10( np.median( civ_t ) )
ovi_t_med = np.log10( np.median( ovi_t ) )


# Metallicity
mg_z_mean = np.log10( np.mean( mgii_z ) )
civ_z_mean = np.log10( np.mean( civ_z ) )
ovi_z_mean = np.log10( np.mean( ovi_z ) )

mg_z_min = np.log10( np.min( mgii_z ) )
civ_z_min = np.log10( np.min( civ_z ) )
ovi_z_min = np.log10( np.min( ovi_z ) )

mg_z_max = np.log10( np.max( mgii_z ) )
civ_z_max = np.log10( np.max( civ_z ) )
ovi_z_max = np.log10( np.max( ovi_z ) )

mg_z_std = np.log10( np.std( mgii_z ) )
civ_z_std = np.log10( np.std( civ_z ) )
ovi_z_std = np.log10( np.std( ovi_z ) )

mg_z_med = np.log10( np.median( mgii_z ) )
civ_z_med = np.log10( np.median( civ_z ) )
ovi_z_med = np.log10( np.median( ovi_z ) )



# Cell size
mg_l_mean = np.mean( mgii_l ) 
civ_l_mean = np.mean( civ_l ) 
ovi_l_mean =  np.mean( ovi_l ) 

mg_l_min = np.min( mgii_l ) 
civ_l_min = np.min( civ_l ) 
ovi_l_min = np.min( ovi_l ) 

mg_l_max =  np.max( mgii_l ) 
civ_l_max = np.max( civ_l ) 
ovi_l_max = np.max( ovi_l ) 

mg_l_std = np.std( mgii_l ) 
civ_l_std = np.std( civ_l ) 
ovi_l_std = np.std( ovi_l ) 

mg_l_med = np.median( mgii_l ) 
civ_l_med = np.median( civ_l ) 
ovi_l_med = np.median( ovi_l ) 





s = 'Cells that do not have HI absorption\n\n'
fo.write(s)

s = 'Typical MgII cell properties\n'
fo.write(s)
s = 'Property \t Min \t Max \t Mean \t Med \t Std Dev \n'
fo.write(s)
s = 'Density \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(mg_nH_min, mg_nH_max, mg_nH_mean, mg_nH_med, mg_nH_std)
fo.write(s)
s = 'Temperature \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(mg_t_min, mg_t_max, mg_t_mean, mg_t_med, mg_t_std)
fo.write(s)
s = 'Metalicity \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(mg_z_min, mg_z_max, mg_z_mean, mg_z_med, mg_z_std)
fo.write(s)
s = 'Cell Size \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(mg_l_min, mg_l_max, mg_l_mean, mg_l_med, mg_l_std)
fo.write(s)
s = 'Fraction of MgII absorbing cells that do not have HI absorption: {0:.2%}\n'.format(mgFrac)
fo.write(s)


s = '\nTypical CIV cell properties\n'
fo.write(s)
s = 'Property \t Min \t Max \t Mean \t Med \t Std Dev \n'
fo.write(s)
s = 'Density \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(civ_nH_min, civ_nH_max, civ_nH_mean, civ_nH_med, civ_nH_std)
fo.write(s)
s = 'Temperature \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(civ_t_min, civ_t_max, civ_t_mean, civ_t_med, civ_t_std)
fo.write(s)
s = 'Metalicity \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(civ_z_min, civ_z_max, civ_z_mean, civ_z_med, civ_z_std)
fo.write(s)
s = 'Cell Size \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(civ_l_min, civ_l_max, civ_l_mean, civ_l_med, civ_l_std)
fo.write(s)
s = 'Fraction of CIV absorbing cells that do not have HI absorption: {0:.2%}\n'.format(civFrac)
fo.write(s)

s = '\nTypical OVI cell properties\n'
fo.write(s)
s = 'Property \t Min \t Max \t Mean \t Med \t Std Dev \n'
fo.write(s)
s = 'Density \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(ovi_nH_min, ovi_nH_max, ovi_nH_mean, ovi_nH_med, ovi_nH_std)
fo.write(s)
s = 'Temperature \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(ovi_t_min, ovi_t_max, ovi_t_mean, ovi_t_med, ovi_t_std)
fo.write(s)
s = 'Metalicity \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(ovi_z_min, ovi_z_max, ovi_z_mean, ovi_z_med, ovi_z_std)
fo.write(s)
s = 'Cell Size \t {0:.2f} \t {1:.2f} \t {2:.2f} \t {3:.2f} \t {4:.2f} \n'.format(ovi_l_min, ovi_l_max, ovi_l_mean, ovi_l_med, ovi_l_std)
fo.write(s)
s = 'Fraction of OVI absorbing cells that do not have HI absorption:  {0:.2%}\n'.format(oviFrac)
fo.write(s)




fo.close()
