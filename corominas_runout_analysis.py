#modules and library
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


#functions
def compute_HoL(V,A,B,CL):
    #A,B,CL regression parameters (intercept, slope, 95% confidence limit) from Corominas 1996, Table 1
    #V landslide volume in m^3 
    
    #HoL from regression line
    logHoL=A+B*np.log10(V)

    #HoL from lower bound 95% confidence interval
    logHoL_lb=A+B*np.log10(V)-CL
    
    return logHoL, logHoL_lb

def plot_HoL_vs_V(V,A,B,CL,LT_P):
    #A,B,CL regression parameters (intercept, slope, 95% confidence limit) from Corominas 1996, Table 1
    #V landslide volume in m^3
    #LT_P landslide type and path

    ax=plt.gca()

    #plotting lines    
    V0=1,max(V)*10
    logHoL,logHoL_lb=compute_HoL(V0,A,B,CL)
    ax.plot(V0,logHoL,'b-',label='mean')
    ax.plot(V0,logHoL_lb,'b:',label='95%,lower')

    #plotting points
    logHoL,logHoL_lb=compute_HoL(V,A,B,CL)
    ax.plot(V,logHoL,'ro')
    ax.plot(V,logHoL_lb,'rs')

    #formatting figure
    ax.semilogx()
    ax.set_xlabel('log V')
    ax.set_ylabel('log H/L')
    ax.set_title(LT_P)
    plt.grid(b=None, which='major',axis='both')

    return logHoL,logHoL_lb

def user_input():#USER INPUTS
    while 1==1:
        try:
            s = raw_input("Input landslide type-path  (separate with comma): ")
            ltp_id = map(int, s.split(','))
            all(isinstance(a, int) for a in ltp_id)
            all(0 <= a < len(df) for a in ltp_id)
            break
        except:
            continue
        
    while 1==1:
        try:
            s = raw_input("Input values of mean thickness of landslide (separate with comma): ")
            thickness = map(float, s.split(','))
            all(isinstance(a, float) for a in thickness)
            try:
                len(thickness)
                thickness=np.array(thickness)
            except:
                break
            break
        except:
            continue

    while 1==1:
        try:
            s = raw_input("Input landslide area: ")
            area = float(s)
            break
        except:
            continue

    while 1==1:
        try:
            s = raw_input("Height of source (scarp) from the base (exposure), in meters: ")
            Elev = float(s)
            break
        except:
            continue

    while 1==1:
        try:
            s = raw_input("Horizontal distance from scarp to exposure, in meters: ")
            Base = float(s)
            break
        except:
            continue

    return ltp_id, thickness , area, Elev, Base   
    
    

    
    

df=pd.read_csv('Corominas1996_table1.csv',header=0)
print df
ltp_id, thickness , area, Elev, Base=user_input()
volume=thickness*area


fig1,ax1=plt.subplots(ncols=len(ltp_id),sharex=True,sharey=True,figsize=(len(ltp_id*3),4))


ax1_in=0
for i in ltp_id:   
    
    Ls,Path,A,B,CL=df.Landslide_type.values[i],df.Path.values[i],df.A.values[i],df.B.values[i],df.CL.values[i]

    #plotting H/L vs V
    plt.sca(ax1[ax1_in])
    logHoL,logHoL_lb=plot_HoL_vs_V(volume,A,B,CL,Ls+' -\n'+Path)
    if ax1_in==len(ax1)-1:plt.gca().legend(fontsize='small')
    ax1_in=ax1_in+1

    
    #plotting runout profile for current landslide type - path
    #logHoL,logHoL_lb=compute_HoL(volume,A,B,CL)
    print i
    print A,B,CL
    print Elev/10**logHoL
    print Elev/10**logHoL_lb

    #plotting runout profile for current landslide type - path

    fig2,ax2=plt.subplots(nrows=len(volume),sharex=True,sharey=True)

    for v_i in range(len(volume)):
        plt.sca(ax2[v_i])
        ax=plt.gca()

       


        #plotting runout profiles (95%, upper)
        H=10**logHoL_lb[v_i]*np.array([0,Base])+Elev
        maxL=Elev/10**logHoL_lb[v_i]
        ax.add_patch(patches.Rectangle((0, 0), maxL, 0.1*Elev, fc='r',lw=0,alpha=0.3,label=('maximum extent (95%='+str(int(round(maxL)))+' $m$)')))
        ax.plot([0,maxL],[Elev,0], 'b:')
        

        #plotting runout profiles (mean)
        H=10**logHoL[v_i]*np.array([0,Base])+Elev
        meanL=Elev/10**logHoL[v_i]
        ax.add_patch(patches.Rectangle((0, 0), meanL, 0.2*Elev, fc='r',lw=0,label=('most likely extent (mean='+str(int(round(meanL)))+' $m$)')))
        plt.plot([0,meanL],[Elev,0], 'b-')
       
        #plotting scarp-exposure profile
        ax.plot([0,Base],[Elev,0],'k--',lw=3)
        ax.plot([0,0,maxL],[Elev,0,0],'k-')
        ax.plot([Base],[0.1*Elev],'y^',markersize=10,label='nearest exposure ('+str(Base)+' $m$)')

        if v_i==len(volume)-1:
            ax.set_xlabel('distance, $m$')
        ax.legend(fontsize='x-small',numpoints=1)
        #ax.yaxis.set_label_position("right")
        ax.set_title('v='+str(int(round(volume[v_i],0)))+' $m^3$, t='+str(thickness[v_i])+' $m$',fontsize='medium')
        plt.setp(ax.get_xticklabels(), fontsize='small')
        plt.setp(ax.get_yticklabels(), fontsize='small')

    ax.axis('equal')
    
    fig2.tight_layout()
    plt.subplots_adjust(top=0.9)
    fig2.suptitle(Ls+' - '+Path)
    
   
fig1.tight_layout()

plt.show()
        
        
    
    


 
