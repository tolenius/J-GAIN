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
        
cases = {
"case1":{
         "wildcard":"6_01_*.vs.6_01_1010.interp-logj2.jump1"
        },
        
"case2":{
         "wildcard":"6_02_*.vs.6_02_1010.interp-logj2.jump1"
        }        
}

for casekey in cases.keys():
    print("case: " + casekey)
    case = cases[casekey]
    file_wildcard = case["wildcard"]
    
    
    root='/nobackup/smhid17/users/sm_danya/MAINLWD/AeroSourses/Article/cases/genpv3-6_*/**/'
    wildcard = os.path.join(root,file_wildcard)

    files = glb.glob(wildcard)
    files.sort()
    plot_ref=False
    print("files count: " + str(len(files)))

    for afile in files[0:2]:
        
        print("working on file: ",afile)
        figure(figsize=(8, 6), dpi=300)
        plt.rcParams['font.size'] = '14'
        ax = plt.gca()
        
        vsindex=afile.find('.interp')
        si=vsindex-2
        ei=vsindex
        ref_pow=afile[si:ei]
        
        
        
        vsindex=afile.find('.vs.')
        si=vsindex-2
        ei=vsindex
        exp_pow=afile[si:ei]
        
        output_file="RE_" + file_wildcard.replace("*", exp_pow+exp_pow)+".png"
        
        ref = "$2^{" + str(int(ref_pow)) + "}+1$"
        exp = "exp_$2^{" + str(int(exp_pow)) + "}+1$"
        diff = "$2^{" + str(int(exp_pow)) + "}+1$"
        df = pd.read_csv(afile)[['H2SO4','NH3','ref','exp']].rename(columns={"ref": ref , "exp": exp  }) 
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
        ax.set_xlabel('NH3 [cm$^{-3}$]')
        ax.set_ylabel('H2SO4 [cm$^{-3}$]')
        
        
        plt.title("Nucleation Relative Error(%).\nNo. of Points: " + "$(2^{" + str(int(exp_pow)) + "}+1)$ x $(2^{" + str(int(exp_pow)) + "}+1)$")
        plt.savefig(output_file, dpi=300)
        
        plt.close()
    
    













