import pandas as pd
import glob2 as glb 
import os 
import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.pyplot import figure
import math 
import seaborn as sns
import math

import matplotlib as mpl
mpl.use('Agg')

colors = [ 
'#4040ff', '#6060ff', '#8080ff', '#9f9fff', '#bfbfff', '#dfdfff', 
'#ffffff', 
'#ffffff', 
'#ffdfdf', '#ffbfbf', '#ff9f9f', '#ff8080', '#ff6060', '#ff4040', '#ff2020', '#dd0000','#aa0000'
]
        
 
ifile="./t2t/t2t.interp-logj2.jump1"  

print("Plotting  on file: ",ifile)
figure(figsize=(8, 6), dpi=300)
plt.rcParams['font.size'] = '14'
ax = plt.gca()


ref_pow=6
exp_pow=4

ref = "$2^{" + str(int(ref_pow)) + "}+1$"
exp = "exp_$2^{" + str(int(exp_pow)) + "}+1$"
diff = "$2^{" + str(int(exp_pow)) + "}+1$"
df = pd.read_csv(ifile)[['H2SO4','NH3','ref','exp']].rename(columns={"ref": ref , "exp": exp  }) 
df['H2SO4'] = df.apply( lambda x : x['H2SO4']/1000000 , axis = 1)
df['NH3'] = df.apply( lambda x : x['NH3']/1000000 , axis = 1)
df['DIFF'] = df.apply( lambda x :100*(x[exp]-x[ref])/x[ref] if x[ref] > 1 and (abs(x[exp]-x[ref]) > 1) else 0 , axis = 1)

 

nn = len(df)

n = int(math.sqrt(nn))
nh = np.array(df['NH3'])
h2 = np.array(df['H2SO4'])
diff = np.array(df['DIFF'])
nh.shape = (n,n)
h2.shape = (n,n)
diff.shape = (n,n)

#sns.kdeplot(data =  , x='H2SO4', y='NH3', z = 'DIFF',  log_scale=(10,10), shade=True, bw_adjust=.5)
plt.xscale('log')
plt.yscale('log')  
cbarticks = np.arange(-7,11,1)        
#im = plt.contourf(nh, h2, diff, cbarticks, alpha=0.7, cmap=plt.cm.jet, vmin = -6, vmax = 10, colors=colors)
im = plt.contourf(nh, h2, diff, cbarticks, alpha=0.7,  vmin = -7, vmax = 10, colors=colors)
cbar = plt.colorbar(im, ax=ax,ticks=cbarticks)
cbar.set_label('Relative Error (%)')
ax.set_xlabel('[NH${_3}$][cm$^{-3}$]')
ax.set_ylabel('[H${_2}$SO$_4$][cm$^{-3}$]')

num=2**int(exp_pow)+1
plt.title("Nucleation Relative Error(%).\nNo. of Points: " + "$(2^{" + str(int(exp_pow)) + "}+1)$x$(2^{" + str(int(exp_pow)) + "}+1) = " + str(num) + " $x$ " + str(num) + "$")
plt.savefig("ref.vs.exp.png", dpi=300)
#plt.show()
plt.close()
    
    













