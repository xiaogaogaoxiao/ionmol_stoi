import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import interpolate
from matplotlib.legend_handler import HandlerLine2D
plt.rc("text", usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{wasysym} \usepackage{stmaryrd}')

convMb=0.28002*100
convcm2=0.28002
co_p='tab:red'
co_ap='tab:red'
co_he='tab:green'
co_be='tab:blue'
co_c='tab:purple'
co_o='tab:orange'

def convkeV(x):
    return 25.0*(x**2)

def convMeV(x):
    return 25.0*(x**2)/1000.

def interp_XS(XY,ix,iy,kind):
    x=XY[ix].values
    y=XY[iy].values
    n=len(x)
    spline = interpolate.interp1d(x, y, kind=kind)
    xnew = np.arange(x[0], x[n-1], 0.01)
    ynew = spline(xnew)
    XYnew = pd.DataFrame({ix:xnew, iy:ynew})
    return XYnew

class HandlerXoffset(HandlerLine2D):
    def __init__(self, marker_pad=0.3, numpoints=1, x_offset=0,  **kw):
        HandlerLine2D.__init__(self, marker_pad=marker_pad, numpoints=numpoints, **kw)
        self._xoffset = x_offset
    def get_xdata(self, legend, xdescent, ydescent, width, height, fontsize):
        numpoints = self.get_numpoints(legend)

        if numpoints > 1:
            # we put some pad here to compensate the size of the
            # marker
            xdata = np.linspace(-xdescent + self._marker_pad * fontsize,
                                width - self._marker_pad * fontsize,
                                numpoints) - self._xoffset
            xdata_marker = xdata
        elif numpoints == 1:
            xdata = np.linspace(-xdescent, width, 2) - self._xoffset
            xdata_marker = [0.5 * width - 0.5 * xdescent - self._xoffset]
        return xdata, xdata_marker

def erp_avg(erp):
    return sum(erp)/len(erp)

def get_error(xexp,yexp,dfint,chmol,ionch):
    xdef=convMeV(dfint.loc[:]['ENERGY'])/ionch**(2-alpha)
    idx=dfint.index[abs(xdef-xexp)<=0.002]
    if len(idx)==0: idx=dfint.index[abs(xdef-xexp)<=0.0059]
    xtheo=convMeV(dfint.iloc[idx]['ENERGY']).tolist()[0]/ionch**(2-alpha)
    ytheo=dfint.iloc[idx][chmol].tolist()[0]*convcm2/ionch**(alpha)
#     erp=abs(yexp-ytheo)/(0.5*(yexp+ytheo))*100
    erp=(yexp-ytheo)/ytheo*100
#     print("   xexp=",xexp,"xtheo=",xtheo)
#     print("   yexp=",yexp,"ytheo=",ytheo)
    return xtheo,ytheo,erp

def interp_XS_EP(xnew,X,Y,kind):
    x=X.values
    y=Y.values
    spline = interpolate.interp1d(x, y, kind=kind, fill_value='extrapolate')
    ynew = spline(xnew)
    return ynew

folder="~/ionmol_stoi/calculos/"

PS1=pd.read_csv(folder+"P/P_CDW_adn_NV.dat",sep=' ',delimiter='\s+',header='infer')
PS2=pd.read_csv(folder+"P/P_CDW_CHn_NV.dat",sep=' ',delimiter='\s+',header='infer')
PS3=pd.read_csv(folder+"P/P_CDW_PYR_NV.dat",sep=' ',delimiter='\s+',header='infer')
HeS1=pd.read_csv(folder+"He2/He2_CDW_adn_NV.dat",sep=' ',delimiter='\s+',header='infer')
HeS2=pd.read_csv(folder+"He2/He2_CDW_CHn_NV.dat",sep=' ',delimiter='\s+',header='infer')
HeS3=pd.read_csv(folder+"He2/He2_CDW_PYR_NV.dat",sep=' ',delimiter='\s+',header='infer')
BeS1=pd.read_csv(folder+"Be4/Be4_CDW_adn_NV.dat",sep=' ',delimiter='\s+',header='infer')
BeS2=pd.read_csv(folder+"Be4/Be4_CDW_CHn_NV.dat",sep=' ',delimiter='\s+',header='infer')
BeS3=pd.read_csv(folder+"Be4/Be4_CDW_PYR_NV.dat",sep=' ',delimiter='\s+',header='infer')
CS1=pd.read_csv(folder+"C6/C6_CDW_adn_NV.dat",sep=' ',delimiter='\s+',header='infer')
CS2=pd.read_csv(folder+"C6/C6_CDW_CHn_NV.dat",sep=' ',delimiter='\s+',header='infer')
CS3=pd.read_csv(folder+"C6/C6_CDW_PYR_NV.dat",sep=' ',delimiter='\s+',header='infer')
OS1=pd.read_csv(folder+"O8/O8_CDW_adn_NV.dat",sep=' ',delimiter='\s+',header='infer')
OS2=pd.read_csv(folder+"O8/O8_CDW_CHn_NV.dat",sep=' ',delimiter='\s+',header='infer')
OS3=pd.read_csv(folder+"O8/O8_CDW_PYR_NV.dat",sep=' ',delimiter='\s+',header='infer')

PNS1=pd.read_csv(folder+"P/P_CDW_adn.dat",sep=' ',delimiter='\s+',header='infer')
PNS2=pd.read_csv(folder+"P/P_CDW_CHn.dat",sep=' ',delimiter='\s+',header='infer')
PNS3=pd.read_csv(folder+"P/P_CDW_PYR.dat",sep=' ',delimiter='\s+',header='infer')
HeNS1=pd.read_csv(folder+"He2/He2_CDW_adn.dat",sep=' ',delimiter='\s+',header='infer')
HeNS2=pd.read_csv(folder+"He2/He2_CDW_CHn.dat",sep=' ',delimiter='\s+',header='infer')
HeNS3=pd.read_csv(folder+"He2/He2_CDW_PYR.dat",sep=' ',delimiter='\s+',header='infer')
BeNS1=pd.read_csv(folder+"Be4/Be4_CDW_adn.dat",sep=' ',delimiter='\s+',header='infer')
BeNS2=pd.read_csv(folder+"Be4/Be4_CDW_CHn.dat",sep=' ',delimiter='\s+',header='infer')
BeNS3=pd.read_csv(folder+"Be4/Be4_CDW_PYR.dat",sep=' ',delimiter='\s+',header='infer')
CNS1=pd.read_csv(folder+"C6/C6_CDW_adn.dat",sep=' ',delimiter='\s+',header='infer')
CNS2=pd.read_csv(folder+"C6/C6_CDW_CHn.dat",sep=' ',delimiter='\s+',header='infer')
CNS3=pd.read_csv(folder+"C6/C6_CDW_PYR.dat",sep=' ',delimiter='\s+',header='infer')
ONS1=pd.read_csv(folder+"O8/O8_CDW_adn.dat",sep=' ',delimiter='\s+',header='infer')
ONS2=pd.read_csv(folder+"O8/O8_CDW_CHn.dat",sep=' ',delimiter='\s+',header='infer')
ONS3=pd.read_csv(folder+"O8/O8_CDW_PYR.dat",sep=' ',delimiter='\s+',header='infer')

head_adn=list(PS1.columns.values)
head_CHn=list(PS2.columns.values)
head_PYR=list(PS3.columns.values)

headNS_adn=list(PNS1.columns.values)
headNS_CHn=list(PNS2.columns.values)
headNS_PYR=list(PNS3.columns.values)

# water data 
# E -- MeV/amu
# CS -- a.u.
PS4=pd.read_csv(folder+"P/P_CDW_H2O_n.txt",sep=' ',delimiter='\s+',header='infer')
HeS4=pd.read_csv(folder+"He2/He2_CDW_H2O_n.txt",sep=' ',delimiter='\s+',header='infer')
BeS4=pd.read_csv(folder+"Be4/Be4_CDW_H2O_n.txt",sep=' ',delimiter='\s+',header='infer')
CS4=pd.read_csv(folder+"C6/C6_CDW_H2O_n.txt",sep=' ',delimiter='\s+',header='infer')
OS4=pd.read_csv(folder+"O8/O8_CDW_H2O_n.txt",sep=' ',delimiter='\s+',header='infer')

PS4['Energy']=np.sqrt(PS4['Energy']*1000/25.)
HeS4['Energy']=np.sqrt(HeS4['Energy']*1000/25.)
BeS4['Energy']=np.sqrt(BeS4['Energy']*1000/25.)
CS4['Energy']=np.sqrt(CS4['Energy']*1000/25.)
OS4['Energy']=np.sqrt(OS4['Energy']*1000/25.)

# Convert energy from keV to a.u.

PS1['ENERGY']=np.sqrt(PS1['ENERGY']/25.)
PS2['ENERGY']=np.sqrt(PS2['ENERGY']/25.)
PS3['ENERGY']=np.sqrt(PS3['ENERGY']/25.)
HeS1['ENERGY']=np.sqrt(HeS1['ENERGY']/25.)
HeS2['ENERGY']=np.sqrt(HeS2['ENERGY']/25.)
HeS3['ENERGY']=np.sqrt(HeS3['ENERGY']/25.)
BeS1['ENERGY']=np.sqrt(BeS1['ENERGY']/25.)
BeS2['ENERGY']=np.sqrt(BeS2['ENERGY']/25.)
BeS3['ENERGY']=np.sqrt(BeS3['ENERGY']/25.)
CS1['ENERGY']=np.sqrt(CS1['ENERGY']/25.)
CS2['ENERGY']=np.sqrt(CS2['ENERGY']/25.)
CS3['ENERGY']=np.sqrt(CS3['ENERGY']/25.)
OS1['ENERGY']=np.sqrt(OS1['ENERGY']/25.)
OS2['ENERGY']=np.sqrt(OS2['ENERGY']/25.)
OS3['ENERGY']=np.sqrt(OS3['ENERGY']/25.)
PNS1['ENERGY']=np.sqrt(PNS1['ENERGY']/25.)
PNS2['ENERGY']=np.sqrt(PNS2['ENERGY']/25.)
PNS3['ENERGY']=np.sqrt(PNS3['ENERGY']/25.)
HeNS1['ENERGY']=np.sqrt(HeNS1['ENERGY']/25.)
HeNS2['ENERGY']=np.sqrt(HeNS2['ENERGY']/25.)
HeNS3['ENERGY']=np.sqrt(HeNS3['ENERGY']/25.)
BeNS1['ENERGY']=np.sqrt(BeNS1['ENERGY']/25.)
BeNS2['ENERGY']=np.sqrt(BeNS2['ENERGY']/25.)
BeNS3['ENERGY']=np.sqrt(BeNS3['ENERGY']/25.)
CNS1['ENERGY']=np.sqrt(CNS1['ENERGY']/25.)
CNS2['ENERGY']=np.sqrt(CNS2['ENERGY']/25.)
CNS3['ENERGY']=np.sqrt(CNS3['ENERGY']/25.)
ONS1['ENERGY']=np.sqrt(ONS1['ENERGY']/25.)
ONS2['ENERGY']=np.sqrt(ONS2['ENERGY']/25.)
ONS3['ENERGY']=np.sqrt(ONS3['ENERGY']/25.)

# Convert cross section per electron to molecular

ne_adn=[36,45,42,37,49,54.5,151]
npts=len(head_adn)
for i in range(1,npts):
    name=head_adn[i]
    PS1[name]=PS1[name]*ne_adn[i-1]
    HeS1[name]=HeS1[name]*ne_adn[i-1]
    BeS1[name]=BeS1[name]*ne_adn[i-1]
    CS1[name]=CS1[name]*ne_adn[i-1]
    OS1[name]=OS1[name]*ne_adn[i-1]
    
ne_CHn=[8,10,12,2*4+6,6*4+6]
npts=len(head_CHn)
for i in range(1,npts):
    name=head_CHn[i]
    PS2[name]=PS2[name]*ne_CHn[i-1]
    HeS2[name]=HeS2[name]*ne_CHn[i-1]
    BeS2[name]=BeS2[name]*ne_CHn[i-1]
    CS2[name]=CS2[name]*ne_CHn[i-1]
    OS2[name]=OS2[name]*ne_CHn[i-1]
    
ne_PYR=[4*4+4+4*2,2*4+7+4,4+5+4,4*5+5+4,4*4+8+4]
npts=len(head_PYR)
for i in range(1,npts):
    name=head_PYR[i]
    PS3[name]=PS3[name]*ne_PYR[i-1]
    HeS3[name]=HeS3[name]*ne_PYR[i-1]
    BeS3[name]=BeS3[name]*ne_PYR[i-1]
    CS3[name]=CS3[name]*ne_PYR[i-1]
    OS3[name]=OS3[name]*ne_PYR[i-1]

# Make cubic spline interpolation on data

PS1int=[]
PNS1int=[]
PS1int=interp_XS(PS1,'ENERGY','SLURACILO','cubic')
PNS1int=interp_XS(PNS1,'ENERGY','URACILO','cubic')
for i in head_adn[2:]:
    dum1=interp_XS(PS1,'ENERGY',i,'cubic')
    PS1int[i]=dum1[i]
for i in headNS_adn[2:]:
    dumn1=interp_XS(PNS1,'ENERGY',i,'cubic')
    PNS1int[i]=dumn1[i]

PS2int=[]
PNS2int=[]
PS2int=interp_XS(PS2,'ENERGY','SLCH4','cubic')
PNS2int=interp_XS(PNS2,'ENERGY','CH4','cubic')
for i in head_CHn[2:]:
    dum2=interp_XS(PS2,'ENERGY',i,'cubic')
    PS2int[i]=dum2[i]
for i in headNS_CHn[2:]:
    dumn2=interp_XS(PNS2,'ENERGY',i,'cubic')
    PNS2int[i]=dumn2[i]

PS3int=[]
PNS3int=[]
PS3int=interp_XS(PS3,'ENERGY','SLC4H4N2','cubic')
PNS3int=interp_XS(PNS3,'ENERGY','C4H4N2','cubic')
for i in head_PYR[2:]:
    dum3=interp_XS(PS3,'ENERGY',i,'cubic')
    PS3int[i]=dum3[i]
for i in headNS_PYR[2:]:
    dumn3=interp_XS(PNS3,'ENERGY',i,'cubic')
    PNS3int[i]=dumn3[i]

HeS1int=[]
HeNS1int=[]
HeS1int=interp_XS(HeS1,'ENERGY','SLURACILO','cubic')
HeNS1int=interp_XS(HeNS1,'ENERGY','URACILO','cubic')
for i in head_adn[2:]:
    dum1=interp_XS(HeS1,'ENERGY',i,'cubic')
    HeS1int[i]=dum1[i]
for i in headNS_adn[2:]:
    dumn1=interp_XS(HeNS1,'ENERGY',i,'cubic')
    HeNS1int[i]=dumn1[i]

HeS2int=[]
HeNS2int=[]
HeS2int=interp_XS(HeS2,'ENERGY','SLCH4','cubic')
HeNS2int=interp_XS(HeNS2,'ENERGY','CH4','cubic')
for i in head_CHn[2:]:
    dum2=interp_XS(HeS2,'ENERGY',i,'cubic')
    HeS2int[i]=dum2[i]
for i in headNS_CHn[2:]:
    dumn2=interp_XS(HeNS2,'ENERGY',i,'cubic')
    HeNS2int[i]=dumn2[i]

HeS3int=[]
HeNS3int=[]
HeS3int=interp_XS(HeS3,'ENERGY','SLC4H4N2','cubic')
HeNS3int=interp_XS(HeNS3,'ENERGY','C4H4N2','cubic')
for i in head_PYR[2:]:
    dum3=interp_XS(HeS3,'ENERGY',i,'cubic')
    HeS3int[i]=dum3[i]
for i in headNS_PYR[2:]:
    dumn3=interp_XS(HeNS3,'ENERGY',i,'cubic')
    HeNS3int[i]=dumn3[i]

BeS1int=[]
BeNS1int=[]
BeS1int=interp_XS(BeS1,'ENERGY','SLURACILO','cubic')
BeNS1int=interp_XS(BeNS1,'ENERGY','URACILO','cubic')
for i in head_adn[2:]:
    dum1=interp_XS(BeS1,'ENERGY',i,'cubic')
    BeS1int[i]=dum1[i]
for i in headNS_adn[2:]:
    dumn1=interp_XS(BeNS1,'ENERGY',i,'cubic')
    BeNS1int[i]=dumn1[i]

BeS2int=[]
BeNS2int=[]
BeS2int=interp_XS(BeS2,'ENERGY','SLCH4','cubic')
BeNS2int=interp_XS(BeNS2,'ENERGY','CH4','cubic')
for i in head_CHn[2:]:
    dum2=interp_XS(BeS2,'ENERGY',i,'cubic')
    BeS2int[i]=dum2[i]
for i in headNS_CHn[2:]:
    dumn2=interp_XS(BeNS2,'ENERGY',i,'cubic')
    BeNS2int[i]=dumn2[i]

BeS3int=[]
BeNS3int=[]
BeS3int=interp_XS(BeS3,'ENERGY','SLC4H4N2','cubic')
BeNS3int=interp_XS(BeNS3,'ENERGY','C4H4N2','cubic')
for i in head_PYR[2:]:
    dum3=interp_XS(BeS3,'ENERGY',i,'cubic')
    BeS3int[i]=dum3[i]
for i in headNS_PYR[2:]:
    dumn3=interp_XS(BeNS3,'ENERGY',i,'cubic')
    BeNS3int[i]=dumn3[i]

CS1int=[]
CNS1int=[]
CS1int=interp_XS(CS1,'ENERGY','SLURACILO','cubic')
CNS1int=interp_XS(CNS1,'ENERGY','URACILO','cubic')
for i in head_adn[2:]:
    dum1=interp_XS(CS1,'ENERGY',i,'cubic')
    CS1int[i]=dum1[i]
for i in headNS_adn[2:]:
    dumn1=interp_XS(CNS1,'ENERGY',i,'cubic')
    CNS1int[i]=dumn1[i]

CS2int=[]
CNS2int=[]
CS2int=interp_XS(CS2,'ENERGY','SLCH4','cubic')
CNS2int=interp_XS(CNS2,'ENERGY','CH4','cubic')
for i in head_CHn[2:]:
    dum2=interp_XS(CS2,'ENERGY',i,'cubic')
    CS2int[i]=dum2[i]
for i in headNS_CHn[2:]:
    dumn2=interp_XS(CNS2,'ENERGY',i,'cubic')
    CNS2int[i]=dumn2[i]

CS3int=[]
CNS3int=[]
CS3int=interp_XS(CS3,'ENERGY','SLC4H4N2','cubic')
CNS3int=interp_XS(CNS3,'ENERGY','C4H4N2','cubic')
for i in head_PYR[2:]:
    dum3=interp_XS(CS3,'ENERGY',i,'cubic')
    CS3int[i]=dum3[i]
for i in headNS_PYR[2:]:
    dumn3=interp_XS(CNS3,'ENERGY',i,'cubic')
    CNS3int[i]=dumn3[i]

OS1int=[]
ONS1int=[]
OS1int=interp_XS(OS1,'ENERGY','SLURACILO','cubic')
ONS1int=interp_XS(ONS1,'ENERGY','URACILO','cubic')
for i in head_adn[2:]:
    dum1=interp_XS(OS1,'ENERGY',i,'cubic')
    OS1int[i]=dum1[i]
for i in headNS_adn[2:]:
    dumn1=interp_XS(ONS1,'ENERGY',i,'cubic')
    ONS1int[i]=dumn1[i]

OS2int=[]
ONS2int=[]
OS2int=interp_XS(OS2,'ENERGY','SLCH4','cubic')
ONS2int=interp_XS(ONS2,'ENERGY','CH4','cubic')
for i in head_CHn[2:]:
    dum2=interp_XS(OS2,'ENERGY',i,'cubic')
    OS2int[i]=dum2[i]
for i in headNS_CHn[2:]:
    dumn2=interp_XS(ONS2,'ENERGY',i,'cubic')
    ONS2int[i]=dumn2[i]

OS3int=[]
ONS3int=[]
OS3int=interp_XS(OS3,'ENERGY','SLC4H4N2','cubic')
ONS3int=interp_XS(ONS3,'ENERGY','C4H4N2','cubic')
for i in head_PYR[2:]:
    dum3=interp_XS(OS3,'ENERGY',i,'cubic')
    OS3int[i]=dum3[i]
for i in headNS_PYR[2:]:
    dumn3=interp_XS(ONS3,'ENERGY',i,'cubic')
    ONS3int[i]=dumn3[i]

PS4int=interp_XS(PS4,'Energy','H2O','cubic')
HeS4int=interp_XS(HeS4,'Energy','H2O','cubic')
BeS4int=interp_XS(BeS4,'Energy','H2O','cubic')
CS4int=interp_XS(CS4,'Energy','H2O','cubic')
OS4int=interp_XS(OS4,'Energy','H2O','cubic')

# Input $\bar{p}$ data, change $\sigma$ units and interpolate

APS1=pd.read_csv(folder+"AP/AP_CDW_adn_NV.dat",sep=' ',delimiter='\s+',header='infer')
APS3=pd.read_csv(folder+"AP/AP_CDW_PYR_NV.dat",sep=' ',delimiter='\s+',header='infer')

if APS1['ENERGY'][0]==25: 
    APS1['ENERGY']=np.sqrt(APS1['ENERGY']/25.)
if APS3['ENERGY'][0]==25: 
    APS3['ENERGY']=np.sqrt(APS3['ENERGY']/25.)
    
head_APS1=list(APS1.columns.values)
head_APS3=list(APS3.columns.values)

APS1int=[]
APS1int=interp_XS(APS1,'ENERGY','SLURACILO','cubic')
for i in head_APS1[2:]:
    dum1=interp_XS(APS1,'ENERGY',i,'cubic')
    APS1int[i]=dum1[i]

APS3int=[]
APS3int=interp_XS(APS3,'ENERGY','SLC4H4N2','cubic')
for i in head_APS3[2:]:
    dum1=interp_XS(APS3,'ENERGY',i,'cubic')
    APS3int[i]=dum1[i]

# Input of experimental data

folderexp='~/ionmol_stoi/experimentos/'

A_p=pd.read_csv(folderexp+'adenina+p_iriki2011.dat',sep=' ',delimiter='\s+',header='infer')
A_e=pd.read_csv(folderexp+'adenine+e_rahman2016.dat',sep=' ',delimiter='\s+',header='infer')
C_e=pd.read_csv(folderexp+'cytosine+e_rahman2016.dat',sep=' ',delimiter='\s+',header='infer')
G_e=pd.read_csv(folderexp+'guanine+e_rahman2016.dat',sep=' ',delimiter='\s+',header='infer')
T_e=pd.read_csv(folderexp+'thymine+e_rahman2016.dat',sep=' ',delimiter='\s+',header='infer')
U_p=pd.read_csv(folderexp+'uracilo+p_itoh2013.dat',sep=' ',delimiter='\s+',header='infer')
PYR_p=pd.read_csv(folderexp+'pirimidina+p_wolff2014',sep=' ',delimiter='\s+',header='infer')
bug17=pd.read_csv(folderexp+'e_bug2017.dat',sep=' ',delimiter='\s+',header='infer')
THF_p=pd.read_csv(folderexp+'THF+p_wang2016.dat',sep=' ',delimiter='\s+',header='infer')
THF_e09=pd.read_csv(folderexp+'THF+e_fuss2009.dat',sep=' ',delimiter='\s+',header='infer')
THF_e19=pd.read_csv(folderexp+'THF+e_wolff2019.dat',sep=' ',delimiter='\s+',header='infer')
U_c4=pd.read_csv(folderexp+'uracil+c4_ref14-15.dat',sep=' ',delimiter='\s+',header=None)
U_c6=pd.read_csv(folderexp+'uracil+c6_ref14-15.dat',sep=' ',delimiter='\s+',header=None)
U_f6=pd.read_csv(folderexp+'uracil+f6_ref14-15.dat',sep=' ',delimiter='\s+',header=None)
U_o6=pd.read_csv(folderexp+'uracil+o6_ref14-15.dat',sep=' ',delimiter='\s+',header=None)
U_of8=pd.read_csv(folderexp+'uracil+o8f8_ref14-15.dat',sep=' ',delimiter='\s+',header=None)

A_e['Ep']=A_e['Ee']*1.8375597
C_e['Ep']=C_e['Ee']*1.8375597
G_e['Ep']=G_e['Ee']*1.8375597

iAe=len(A_e[A_e['Ep']>500])
iCe=len(C_e[C_e['Ep']<500])
iGe=len(G_e[G_e['Ep']>500])
iTe=len(T_e[T_e['Ep']>500])
iTHFe19=len(THF_e19[THF_e19['Ep']>500])

# Units 
#   Energy        -- Mev/amu
#   Cross section -- 10^-16 cm^2
H2O_p86=pd.read_csv(folderexp+'water+H1_rudd86.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
H2O_p07=pd.read_csv(folderexp+'water+H1_luna07.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
H2O_he2=pd.read_csv(folderexp+'water+He2_rudd85.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
H2O_li3=pd.read_csv(folderexp+'water+li3_luna2016.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
H2O_c6=pd.read_csv(folderexp+'water+C6_dalcapello09.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
H2O_o8=pd.read_csv(folderexp+'water+O8_tribedi2016.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')

gasesH1=pd.read_csv(folderexp+'gases+H1_rudd83.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
gasesHe2=pd.read_csv(folderexp+'gases+He2_rudd85.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
CH4_he19=pd.read_csv(folderexp+'CH4+He2_luna19.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
CH4_he85=pd.read_csv(folderexp+'CH4+He2_rudd85.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')
CH4_p=pd.read_csv(folderexp+'CH4+H1_rudd83.dat',sep=' ',delimiter='\s+',skiprows=4,header='infer')

# C+6 on adenine 
## Units: 
# - Energy = MeV/amu
# - Cross section = $10^{-16}$ cm$^2$

A_C=[3.5,60.7]

# New scaling:

alpha=1.2

dum=PS1int[convkeV(PS1int['ENERGY']) < 100.]
ie100=len(dum)

dum=PS1int[convkeV(PS1int['ENERGY']) < 50.]
iPe100=len(dum)
dum=HeS1int[convkeV(HeS1int['ENERGY'])/2**(2-alpha) < 50.]
iHe100=len(dum)
dum=BeS1int[convkeV(BeS1int['ENERGY'])/4**(2-alpha) < 50.]
iBe100=len(dum)
dum=CS1int[convkeV(CS1int['ENERGY'])/6**(2-alpha) < 50.]
iCe100=len(dum)
dum=OS1int[convkeV(OS1int['ENERGY'])/8**(2-alpha) < 50.]
iOe100=len(dum)

itickssize=28
imolsize=30
ilabelsize=30
ilegsize=26

iwidth=2
xmin=0.035
xmax=15
ymin=0.5
ymax=3.*45
xtext=0.07
ytext=0.09*45
icharge=0.07
ims=16
imew=2

fig = plt.figure(figsize=(16,16))

### ADENINA ###########################################################
ax1 = plt.subplot(421)
plt.text(xtext, ytext, 'Adenine',fontsize=imolsize)
plt.text(xtext*1.05, ytext*0.5, r'C$_5$H$_5$N$_5$',fontsize=imolsize)

plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/8**(2-alpha),OS1int['SLADENINA'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth,label=r'O$^{+8}$')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/6**(2-alpha),CS1int['SLADENINA'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth,label=r'C$^{+6}$')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/4**(2-alpha),BeS1int['SLADENINA'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth,label=r'Be$^{+4}$')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/2**(2-alpha),HeS1int['SLADENINA'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth,label=r'He$^{+2}$')
plt.plot(convMeV(PS1int['ENERGY'][iPe100:]),PS1int['SLADENINA'][iPe100:]*convcm2,
         co_p,linewidth=iwidth,label=r'H$^{+}$')

legend1=ax1.legend(loc='upper left',bbox_to_anchor=(0.1, 1.075, 0, 0),fontsize=ilegsize,ncol=5)

# experimentos #
Ae, = ax1.plot(A_e['Ep'][0:iAe]/1000.,A_e['CS'][0:iAe],color='dimgray',
               marker='s',markersize=ims,markeredgewidth=imew,markerfacecolor='white',markevery=6,
               linestyle=' ',label='_nolegend_')
Ap, = ax1.plot(A_p['E']/1000.,A_p['CS'],color=co_p,marker='o',
               markersize=ims,markeredgewidth=imew,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
AC, = ax1.plot(A_C[0]/6**(2-alpha),A_C[1]/6**(alpha),color=co_c,marker='d',
                markersize=ims,markeredgewidth=imew,markerfacecolor='white',
                linestyle=' ',label='_nolegend_')

handles, labels = ax1.get_legend_handles_labels();
legend1=ax1.legend(handles[::-1], labels[::-1],loc='upper left',bbox_to_anchor=(0.155,1.32,0,0),fontsize=ilegsize,ncol=6)
ax1.legend([Ap,Ae,AC],['H$^{+}$','e$^{-}$',r'C$^{+6}$'],loc='best',fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False)
ax1.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, left=True, right=True);
ax1.tick_params(direction='in',which='minor',length=5,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
plt.gca().add_artist(legend1)
# ax1.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

### CITOSINA ###########################################################
ax2 = plt.subplot(422)
plt.text(xtext, ytext, 'Cytosine',fontsize=imolsize)
plt.text(xtext*0.95, ytext*0.5, r'C$_4$H$_5$N$_3$O',fontsize=imolsize)

plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/8**(2-alpha),OS1int['SLCITOSINA'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/6**(2-alpha),CS1int['SLCITOSINA'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/4**(2-alpha),BeS1int['SLCITOSINA'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/2**(2-alpha),HeS1int['SLCITOSINA'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS1int['ENERGY'][iPe100:]),PS1int['SLCITOSINA'][iPe100:]*convcm2,co_p,linewidth=iwidth, 
         label='_nolegend_')

# experimentos #
Cite,=plt.plot(C_e['Ep'][iCe:]/1000,C_e['CS'][iCe:],color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew,markerfacecolor='white',markevery=4,
         linestyle=' ',label='e')

ax2.legend([Cite],['e$^{-}$'],loc='best',fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False)
ax2.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, 
                left=True, right=True);
ax2.tick_params(direction='in',which='minor',length=5,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
# ax2.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

### GUANINA ###########################################################
ax3 = plt.subplot(423)
plt.text(xtext, ytext, 'Guanine',fontsize=imolsize)
plt.text(xtext*0.95, ytext*0.5, r'C$_5$H$_5$N$_5$O',fontsize=imolsize)

plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/8**(2-alpha),OS1int['SLGUANINA'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/6**(2-alpha),CS1int['SLGUANINA'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/4**(2-alpha),BeS1int['SLGUANINA'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/2**(2-alpha),HeS1int['SLGUANINA'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS1int['ENERGY'][iPe100:]),PS1int['SLGUANINA'][iPe100:]*convcm2,co_p,linewidth=iwidth,
         label='_nolegend_')

# experimentos #
Guane,=plt.plot(G_e['Ep'][0:iGe]/1000,G_e['CS'][0:iGe],color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew,markerfacecolor='white',markevery=4,
         linestyle=' ')

ax3.legend([Guane],['e$^{-}$'],loc='best',fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False)
ax3.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, 
                left=True, right=True);
ax3.tick_params(direction='in',which='minor',length=5,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
ax3.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))
# ax3.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

### TIMINA ###########################################################
ax4 = plt.subplot(424)
plt.text(xtext, ytext, 'Thymine',fontsize=imolsize)
plt.text(xtext*0.95, ytext*0.5, r'C$_5$H$_6$N$_2$O$_2$',fontsize=imolsize)

plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/8**(2-alpha),OS1int['SLTIMINA'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/6**(2-alpha),CS1int['SLTIMINA'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/4**(2-alpha),BeS1int['SLTIMINA'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/2**(2-alpha),HeS1int['SLTIMINA'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS1int['ENERGY'][iPe100:]),PS1int['SLTIMINA'][iPe100:]*convcm2,co_p,linewidth=iwidth,
         label='_nolegend_')

# experimentos #
Tyme, = plt.plot(T_e['Ep'][0:iTe]/1000,T_e['CS'][0:iTe],color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew,markerfacecolor='white',markevery=4,
         linestyle=' ')

ax4.legend([Tyme],['e$^{-}$'],loc='best',fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False)
ax4.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True,
                left=True, right=True);
ax4.tick_params(direction='in',which='minor',length=5,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
ax4.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

### URACILO ###########################################################
ax5 = plt.subplot(425)
plt.text(xtext, ytext, 'Uracil',fontsize=imolsize)
plt.text(xtext*0.8, ytext*0.5, r'C$_4$H$_4$N$_2$O$_2$',fontsize=imolsize)

plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/8**(2-alpha),OS1int['SLURACILO'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth,label=r'O$^{+8}$')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/6**(2-alpha),CS1int['SLURACILO'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth,label=r'C$^{+6}$')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/4**(2-alpha),BeS1int['SLURACILO'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth,label=r'Be$^{+4}$')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/2**(2-alpha),HeS1int['SLURACILO'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth,label=r'He$^{+2}$')
plt.plot(convMeV(PS1int['ENERGY'][iPe100:]),PS1int['SLURACILO'][iPe100:]*convcm2,co_p,linewidth=iwidth,
         label=r'H$^{+}$')

# experiments
Up, = ax5.plot(U_p['E']/1000,U_p['CS'],color=co_p,
               marker='^',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
               label='_nolegend_')
Uc4, = ax5.plot(U_c4[1]/1000/4**(2-alpha),U_c4[3]*4**(2-alpha),color=co_be,
                marker=r'$\ominus$',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',
                linestyle=' ',label='_nolegend_')
Uc6, = ax5.plot(U_c6[1]/1000/6**(2-alpha),U_c6[2]*6**(2-alpha),color=co_c,
                marker=r'$\oplus$',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',
                linestyle=' ',label='_nolegend_')
Uo6, = ax5.plot(U_o6[1]/1000/6**(2-alpha),U_o6[3]*6**(2-alpha),color=co_c,
                marker=r'$\oplus$',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',
                linestyle=' ',label='_nolegend_')
Uf6, = ax5.plot(U_f6[1]/1000/6**(2-alpha),U_f6[2]*6**(2-alpha),color=co_c,
                marker=r'$\oplus$',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',
                linestyle=' ',label='_nolegend_')
Uof8, = ax5.plot(U_of8[0]/1000/8**(2-alpha),U_of8[1]*8**(2-alpha),color=co_o,
                 marker=r'$\otimes$',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')

ax5.legend([Up,Uc4,Uc6,Uof8],['H$^{+}$','C$^{+4}$','C$^{+6}$, O$^{+6}$, F$^{+6}$','O$^{+6}$, F$^{+6}$'],
           loc='upper right',fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False)

plt.ylabel(r"Ionization Cross Section/$Z^{\alpha}$ ($10^{-16}$ cm$^2$)", fontsize=ilabelsize,
           labelpad=15,position=(1,1))
ax5.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True,
                left=True, right=True);
ax5.tick_params(direction='in',which='minor',length=6,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
# ax5.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

### PYRIMIDINA ###########################################################
ax6 = plt.subplot(426)
plt.text(xtext, ytext, 'Pyrimidine',fontsize=imolsize)
plt.text(xtext*1.15, ytext*0.5, r'C$_4$H$_4$N$_2$',fontsize=imolsize)

plt.plot(convMeV(OS3int['ENERGY'][iOe100:])/8**(2-alpha),OS3int['SLC4H4N2'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS3int['ENERGY'][iCe100:])/6**(2-alpha),CS3int['SLC4H4N2'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS3int['ENERGY'][iBe100:])/4**(2-alpha),BeS3int['SLC4H4N2'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS3int['ENERGY'][iHe100:])/2**(2-alpha),HeS3int['SLC4H4N2'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS3int['ENERGY'][iPe100:]),PS3int['SLC4H4N2'][iPe100:]*convcm2,co_p,linewidth=iwidth,
         label='_nolegend_')

# experiments
PYRe, = ax6.plot(bug17['Ep(keV)']/1000,bug17['e+PY'],color='dimgray',
                 marker='>',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
                 label='p')
PYRp, = ax6.plot(PYR_p['E']/1000,PYR_p['au']*convcm2,color=co_p,
                 marker='v',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
                 label='p')

ax6.legend([PYRp,PYRe],['H$^{+}$','e$^{-}$'],loc='best',fontsize=22,handletextpad=0.01,labelspacing=0.1,
           frameon=False)
ax6.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True,
                left=True, right=True);
ax6.tick_params(direction='in',which='minor',length=6,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')

### THF ###########################################################
ax7 = plt.subplot(427)
plt.text(xtext, ytext, 'THF',fontsize=imolsize)
plt.text(xtext*0.9, ytext*0.5, r'C$_4$H$_8$O',fontsize=imolsize)

plt.plot(convMeV(OS3int['ENERGY'][iOe100:])/8**(2-alpha),OS3int['SLC4H8O1'][iOe100:]*convcm2/8**(alpha),
         co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS3int['ENERGY'][iCe100:])/6**(2-alpha),CS3int['SLC4H8O1'][iCe100:]*convcm2/6**(alpha),
         co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS3int['ENERGY'][iBe100:])/4**(2-alpha),BeS3int['SLC4H8O1'][iBe100:]*convcm2/4**(alpha),
         co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS3int['ENERGY'][iHe100:])/2**(2-alpha),HeS3int['SLC4H8O1'][iHe100:]*convcm2/2**(alpha),
         co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS3int['ENERGY'][iPe100:]),PS3int['SLC4H8O1'][iPe100:]*convcm2,co_p,linewidth=iwidth,
         label='_nolegend_')

# experiments
THFe09, = ax7.plot(THF_e09['Ep(kev)']/1000,THF_e09['CS'],color='dimgray',
                   marker='*',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
                   label='p')
THFe17, = ax7.plot(bug17['Ep(keV)']/1000,bug17['e+THF'],color='dimgray',
                   marker='<',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
                   label='p')
THFe19, = ax7.plot(convMeV(THF_e19['v(au)'][0:iTHFe19]),THF_e19['au'][0:iTHFe19]*convcm2,color='dimgray',
                   marker='>',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
                   label='p')
THFp, = ax7.plot(THF_p['Ep']/1000,THF_p['au']*convcm2,color=co_p,
                 marker='D',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',
                 label='p')

ax7.legend([THFp,THFe09,THFe17,THFe19],['H$^+$','e$^{-}$','e$^{-}$','e$^{-}$'],loc='upper right',
           fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False)
ax7.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, 
                left=True, right=True);
ax7.tick_params(direction='in',which='minor',length=6,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
ax7.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))
# ax7.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

### WATER ###########################################################
ax8 = plt.subplot(428)
plt.text(xtext*0.9, ytext, 'Water',fontsize=imolsize)
plt.text(xtext*1.2, ytext*0.5, r'H$_2$O',fontsize=imolsize)
plt.text(xtext*2.7, ytext*0.75, r'$\times\,5$',fontsize=imolsize)

plt.plot(convMeV(OS4int['Energy'][iOe100:])/8**(2-alpha),OS4int['H2O'][iOe100:]*convcm2/8**(alpha)*5,
         co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS4int['Energy'][iCe100:])/6**(2-alpha),CS4int['H2O'][iCe100:]*convcm2/6**(alpha)*5,
         co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS4int['Energy'][iBe100:])/4**(2-alpha),BeS4int['H2O'][iBe100:]*convcm2/4**(alpha)*5,
         co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS4int['Energy'][iHe100:])/2**(2-alpha),HeS4int['H2O'][iHe100:]*convcm2/2**(alpha)*5,
         co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS4int['Energy'][iPe100:])/1**(2-alpha),PS4int['H2O'][iPe100:]*convcm2/1**(alpha)*5,
         co_p,linewidth=iwidth, label='_nolegend_')

# experiments
H2Op86, = ax8.plot(H2O_p86['E'],H2O_p86['TCS']*5,color=co_p,
                    marker='p',markersize=ims+1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
ax8.plot(H2O_p07['E'],H2O_p07['TCS']*5,color='white',
         marker='o',markersize=ims,markeredgewidth=imew-1,markerfacecolor='white',linestyle=' ',label='__none__')
H2Op07, = ax8.plot(H2O_p07['E'],H2O_p07['TCS']*5,color=co_p,
                    marker=r'$\varhexagon$',markersize=ims+2,markeredgewidth=imew-1,linestyle=' ',label='__none__')
ax8.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)*5,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Ohe, = ax8.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)*5,color=co_he,
                    marker=r'$\odot$',markersize=ims+1,markeredgewidth=imew,linestyle=' ',label='__none__')
ax8.plot(H2O_li3['E']/3**(2-alpha),H2O_li3['TCS']/3**(alpha)*5,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oli16, = ax8.plot(H2O_li3['E']/3**(2-alpha),H2O_li3['TCS']/3**(alpha)*5,color='tab:olive',
                    marker=r'$\logof$',markersize=ims+1,markeredgewidth=imew,linestyle=' ',label='__none__')
ax8.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)*5,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oc, = ax8.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)*5,color=co_c,
                    marker=r'$\ovee$',markersize=ims+2,markeredgewidth=imew,linestyle=' ',label='__none__')
ax8.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)*5,color='white',
         marker='s',markersize=ims-3,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oo, = ax8.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)*5,color=co_o,
                    marker=r'$\varobslash$',markersize=ims+1,markeredgewidth=imew,linestyle=' ',label='__none__')

ax8.legend([(H2Op86,H2Op07),H2Ohe,H2Oli16,H2Oc,H2Oo],[r'H$^+$',r'He$^{+2}$',r'Li$^{+3}$','C$^{+6}$','O$^{+8}$'],
           loc='upper right',fontsize=22,handletextpad=0.01,labelspacing=0.1,frameon=False,
           handler_map={H2Op86:HandlerXoffset(x_offset=22,numpoints=1),H2Op07:HandlerXoffset(x_offset=0,numpoints=1)})
plt.xlabel(r"Impact Energy/$Z^{2-\alpha}$ (MeV/amu)", fontsize=ilabelsize,labelpad=10,position=(0,0))
ax8.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, left=True, right=True);
ax8.tick_params(direction='in',which='minor',length=5,bottom=True, top=True, left=True, right=True);

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.yscale('log')
plt.xscale('log')
ax8.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

for ax in fig.get_axes():
    ax.label_outer()
fig.subplots_adjust(hspace=0.0, wspace=0.0)

plt.subplots_adjust(left=0.1, right=0.975, top=0.925, bottom=0.075)
plt.savefig('../zscale_alpha.eps')
#plt.show()


# Universal curve

Zch=[1,2,4,6,8]

ymin=0.02
ymax=2
xmin=0.02
xmax=15
xtext=5
ytext=10
icharge=65/1000
iwidth=0.5
ims=18
imew=0.9

fig = plt.figure(figsize=(16,14))

ax1 = plt.subplot(111)
plt.text(5., 1.2, r'$\alpha=1.2$',fontsize=imolsize+4)
plt.plot(convMeV(PS1int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS1int['SLURACILO'][iPe100:]*convcm2/Zch[0]**alpha/ne_adn[0],
         co_p,linewidth=iwidth, label=r'H$^+$')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS1int['SLURACILO'][iHe100:]*convcm2/Zch[1]**alpha/ne_adn[0],
         co_he,linewidth=iwidth, label=r'He$^{+2}$')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS1int['SLURACILO'][iBe100:]*convcm2/Zch[2]**alpha/ne_adn[0],
         co_be,linewidth=iwidth, label=r'Be$^{+4}$')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS1int['SLURACILO'][iCe100:]*convcm2/Zch[3]**alpha/ne_adn[0],
         co_c,linewidth=iwidth, label=r'C$^{+6}$')
plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS1int['SLURACILO'][iOe100:]*convcm2/Zch[4]**alpha/ne_adn[0],
         co_o,linewidth=iwidth, label=r'O$^{+8}$')

ax1.legend(loc='upper left',bbox_to_anchor=(0.09, 1.09, 0, 0),fontsize=ilegsize,ncol=5)

it=0
for i in head_adn[1:]:
    plt.plot(convMeV(PS1int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS1int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_adn[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS1int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_adn[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS1int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_adn[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS1int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_adn[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS1int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_adn[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

it=0
for i in head_CHn[1:]:
    plt.plot(convMeV(PS2int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS2int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_CHn[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS2int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS2int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_CHn[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS2int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS2int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_CHn[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS2int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS2int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_CHn[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS2int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS2int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_CHn[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

it=0
for i in head_PYR[1:]:
    plt.plot(convMeV(PS3int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS3int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_PYR[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS3int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS3int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_PYR[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS3int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS3int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_PYR[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS3int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS3int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_PYR[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS3int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS3int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_PYR[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

# water
plt.plot(convMeV(OS4int['Energy'][iOe100:])/8**(2-alpha),OS4int['H2O'][iOe100:]*convcm2/8**(alpha)/6,co_o,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS4int['Energy'][iCe100:])/6**(2-alpha),CS4int['H2O'][iCe100:]*convcm2/6**(alpha)/6,co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS4int['Energy'][iBe100:])/4**(2-alpha),BeS4int['H2O'][iBe100:]*convcm2/4**(alpha)/6,co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS4int['Energy'][iHe100:])/2**(2-alpha),HeS4int['H2O'][iHe100:]*convcm2/2**(alpha)/6,co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(PS4int['Energy'][iPe100:])/1**(2-alpha),PS4int['H2O'][iPe100:]*convcm2/1**(alpha)/6,co_p,linewidth=iwidth, label='_nolegend_')

# experiments:

# electrons
Ae, = plt.plot(A_e['Ep'][0:iAe]/1000.,A_e['CS'][0:iAe]/45,color='dimgray',
               marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=6,
               linestyle=' ',label='_nolegend_')
Cite,=plt.plot(C_e['Ep'][iCe:]/1000,C_e['CS'][iCe:]/37,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ',label='e')
Tyme, = plt.plot(T_e['Ep'][0:iTe]/1000,T_e['CS'][0:iTe]/42,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ')
Guane,=plt.plot(G_e['Ep'][0:iGe]/1000,G_e['CS'][0:iGe]/49,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ')
PYRe, = plt.plot(bug17['Ep(keV)']/1000,bug17['e+PY']/28,color='dimgray',
                 marker='>',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
THFe09, = plt.plot(THF_e09['Ep(kev)']/1000,THF_e09['CS']/28,color='dimgray',
                   marker='*',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')
THFe17, = plt.plot(bug17['Ep(keV)']/1000,bug17['e+THF']/28,color='dimgray',
                   marker='<',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')
THFe19, = plt.plot(convMeV(THF_e19['v(au)'][0:iTHFe19]),THF_e19['au'][0:iTHFe19]/28*convcm2,color='dimgray',
                   marker='>',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')

# # gases + H1+
plt.plot(gasesH1['E'],gasesH1['N2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
N2p, = plt.plot(gasesH1['E'],gasesH1['N2'],color=co_p,
                marker=r'$\boxast$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['O2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
O2p, = plt.plot(gasesH1['E'],gasesH1['O2'],color=co_p,
                marker=r'$\boxbox$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['CO'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
COp, = plt.plot(gasesH1['E'],gasesH1['CO'],color=co_p,
                marker=r'$\boxbar$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['CO2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CO2p, = plt.plot(gasesH1['E'],gasesH1['CO2'],color=co_p,
                 marker=r'$\boxbslash$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(CH4_p['E'],CH4_p['TCS']/8,color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4p, = plt.plot(CH4_p['E'],CH4_p['TCS']/8,color=co_p,
                   marker=r'$\boxdot$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# gases + He2+
plt.plot(gasesHe2['E'],gasesHe2['N2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
N2he, = plt.plot(gasesHe2['E'],gasesHe2['N2'],color=co_he,
                 marker=r'$\boxast$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['CO'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
COhe, = plt.plot(gasesHe2['E'],gasesHe2['CO'],color=co_he,
                 marker=r'$\boxbar$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['O2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
O2he, = plt.plot(gasesHe2['E'],gasesHe2['O2'],color=co_he,
                 marker=r'$\boxbox$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['CO2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CO2he, = plt.plot(gasesHe2['E'],gasesHe2['CO2'],color=co_he,
                  marker=r'$\boxbslash$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(CH4_he85['E']/2**(2-alpha),CH4_he85['TCS']/8/2**(alpha),color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4he85, = plt.plot(CH4_he85['E']/2**(2-alpha),CH4_he85['TCS']/8/2**(alpha),color=co_he,
                    marker=r'$\boxdot$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# # methane
plt.plot(CH4_he19['E'],CH4_he19['TCS']/8,color=co_he,
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4he19, = plt.plot(CH4_he19['E'],CH4_he19['TCS']/8,color=co_p,
                    marker=r'$\varowedge$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# water
H2Op86, = plt.plot(H2O_p86['E'],H2O_p86['TCS']/6,color=co_p,
                   marker='p',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',linestyle=' ',label='__nolegend__')
H2Op07, = plt.plot(H2O_p07['E'],H2O_p07['TCS']/6,color=co_p,
                   marker='H',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',linestyle=' ',label='__nolegend__')
plt.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)/6,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Ohe, = plt.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)/6,color=co_he,
                  marker=r'$\odot$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')
plt.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)/6,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oc, = plt.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)/6,color=co_c,
                 marker=r'$\ovee$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')
plt.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)/6,color='white',
         marker='s',markersize=ims-3,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oo, = plt.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)/6,color=co_o,
                 marker=r'$\varobslash$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')

# DNA and RNA
Ap, = plt.plot(A_p['E']/1000,A_p['CS']/45,color=co_p,
               marker='o',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
Up, = plt.plot(U_p['E']/1000,U_p['CS']/36,color=co_p,marker='^',
               markersize=ims,markeredgewidth=2,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
PYRp, = plt.plot(PYR_p['E']/1000,PYR_p['au']/28*convcm2,color=co_p,
                 marker='v',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
THFp, = plt.plot(THF_p['Ep']/1000,THF_p['au']/28*convcm2,color=co_p,
                 marker='D',markersize=ims-2,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
AC, = plt.plot(A_C[0]/6**(2-alpha),A_C[1]/6**(alpha)/45,color=co_c,
               marker='d',markersize=ims,markeredgewidth=imew*2.5,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')

plt.xlabel(r"Impact Energy/Z$^{2-\alpha}$ (MeV/amu)", fontsize=ilabelsize,labelpad=5)
plt.ylabel(r"$\sigma_e/Z^{\alpha}$ ($10^{-16}$ cm$^2$)", fontsize=ilabelsize,labelpad=15)
ax1.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, left=True, right=True);
ax1.tick_params(direction='in',which='minor',length=6,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax) 
plt.yscale('log')
plt.xscale('log')
ax1.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))
# ax1.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

plt.subplots_adjust(left=0.1, right=0.975)

plt.savefig("../zmol_alpha.eps")
#plt.show()


scaledPS1int =PS1int.copy()
scaledHeS1int=HeS1int.copy()
scaledBeS1int=BeS1int.copy()
scaledCS1int =CS1int.copy()
scaledOS1int =OS1int.copy()
scaledPS2int =PS2int.copy()
scaledHeS2int=HeS2int.copy()
scaledBeS2int=BeS2int.copy()
scaledCS2int =CS2int.copy()
scaledOS2int =OS2int.copy()
scaledPS3int =PS3int.copy()
scaledHeS3int=HeS3int.copy()
scaledBeS3int=BeS3int.copy()
scaledCS3int =CS3int.copy()
scaledOS3int =OS3int.copy()
scaledPS4int =PS4int.copy()
scaledHeS4int=HeS4int.copy()
scaledBeS4int=BeS4int.copy()
scaledCS4int =CS4int.copy()
scaledOS4int =OS4int.copy()

it=0
for i in head_adn[1:]:
    scaledPS1int[i] =PS1int[i]*convcm2/Zch[0]**alpha/ne_adn[it]
    scaledHeS1int[i]=HeS1int[i]*convcm2/Zch[1]**alpha/ne_adn[it]
    scaledBeS1int[i]=BeS1int[i]*convcm2/Zch[2]**alpha/ne_adn[it]
    scaledCS1int[i] =CS1int[i]*convcm2/Zch[3]**alpha/ne_adn[it]
    scaledOS1int[i] =OS1int[i]*convcm2/Zch[4]**alpha/ne_adn[it]
    it=it+1
it=0
for i in head_CHn[1:]:
    scaledPS2int[i] =PS2int[i]*convcm2/Zch[0]**alpha/ne_CHn[it]
    scaledHeS2int[i]=HeS2int[i]*convcm2/Zch[1]**alpha/ne_CHn[it]
    scaledBeS2int[i]=BeS2int[i]*convcm2/Zch[2]**alpha/ne_CHn[it]
    scaledCS2int[i] =CS2int[i]*convcm2/Zch[3]**alpha/ne_CHn[it]
    scaledOS2int[i] =OS2int[i]*convcm2/Zch[4]**alpha/ne_CHn[it]
    it=it+1

it=0
for i in head_PYR[1:]:
    scaledPS3int[i] =PS3int[i]*convcm2/Zch[0]**alpha/ne_PYR[it]
    scaledHeS3int[i]=HeS3int[i]*convcm2/Zch[1]**alpha/ne_PYR[it]
    scaledBeS3int[i]=BeS3int[i]*convcm2/Zch[2]**alpha/ne_PYR[it]
    scaledCS3int[i] =CS3int[i]*convcm2/Zch[3]**alpha/ne_PYR[it]
    scaledOS3int[i] =OS3int[i]*convcm2/Zch[4]**alpha/ne_PYR[it]
    it=it+1
    

scaledPS1int['SUM'] =scaledPS1int[head_adn[1:]].sum(axis=1)#/(len(head_adn[1:]))
scaledHeS1int['SUM']=scaledHeS1int[head_adn[1:]].sum(axis=1)#/(len(head_adn[1:]))
scaledBeS1int['SUM']=scaledBeS1int[head_adn[1:]].sum(axis=1)#/(len(head_adn[1:]))
scaledCS1int['SUM'] =scaledCS1int[head_adn[1:]].sum(axis=1)#/(len(head_adn[1:]))
scaledOS1int['SUM'] =scaledOS1int[head_adn[1:]].sum(axis=1)#/(len(head_adn[1:]))

scaledPS2int['SUM'] =scaledPS2int[head_CHn[1:]].sum(axis=1)#/(len(head_CHn[1:]))
scaledHeS2int['SUM']=scaledHeS2int[head_CHn[1:]].sum(axis=1)#/(len(head_CHn[1:]))
scaledBeS2int['SUM']=scaledBeS2int[head_CHn[1:]].sum(axis=1)#/(len(head_CHn[1:]))
scaledCS2int['SUM'] =scaledCS2int[head_CHn[1:]].sum(axis=1)#/(len(head_CHn[1:]))
scaledOS2int['SUM'] =scaledOS2int[head_CHn[1:]].sum(axis=1)#/(len(head_CHn[1:]))

scaledPS3int['SUM'] =scaledPS3int[head_PYR[1:]].sum(axis=1)#/(len(head_PYR[1:]))
scaledHeS3int['SUM']=scaledHeS3int[head_PYR[1:]].sum(axis=1)#/(len(head_PYR[1:]))
scaledBeS3int['SUM']=scaledBeS3int[head_PYR[1:]].sum(axis=1)#/(len(head_PYR[1:]))
scaledCS3int['SUM'] =scaledCS3int[head_PYR[1:]].sum(axis=1)#/(len(head_PYR[1:]))
scaledOS3int['SUM'] =scaledOS3int[head_PYR[1:]].sum(axis=1)#/(len(head_PYR[1:]))

scaledPS4int['H2O'] =PS4int['H2O']*convcm2/1**(alpha)/6
scaledHeS4int['H2O']=HeS4int['H2O']*convcm2/2**(alpha)/6
scaledBeS4int['H2O']=BeS4int['H2O']*convcm2/4**(alpha)/6
scaledCS4int['H2O'] =CS4int['H2O']*convcm2/6**(alpha)/6
scaledOS4int['H2O'] =OS4int['H2O']*convcm2/8**(alpha)/6

resP=scaledPS1int.copy()
resP.drop(head_adn[1:],axis=1,inplace=True)
resHe=scaledHeS1int.copy()
resHe.drop(head_adn[1:],axis=1,inplace=True)
resBe=scaledBeS1int.copy()
resBe.drop(head_adn[1:],axis=1,inplace=True)
resC=scaledCS1int.copy()
resC.drop(head_adn[1:],axis=1,inplace=True)
resO=scaledOS1int.copy()
resO.drop(head_adn[1:],axis=1,inplace=True)
resP['SUM']=scaledPS1int['SUM']+scaledPS2int['SUM']+scaledPS3int['SUM']+scaledPS4int['H2O']
resHe['SUM']=scaledHeS1int['SUM']+scaledHeS2int['SUM']+scaledHeS3int['SUM']+scaledHeS4int['H2O']
resBe['SUM']=scaledBeS1int['SUM']+scaledBeS2int['SUM']+scaledBeS3int['SUM']+scaledBeS4int['H2O']
resC['SUM']=scaledCS1int['SUM']+scaledCS2int['SUM']+scaledCS3int['SUM']+scaledCS4int['H2O']
resO['SUM']=scaledOS1int['SUM']+scaledOS2int['SUM']+scaledOS3int['SUM']+scaledOS4int['H2O']

ntargets=len(head_adn[1:])+len(head_CHn[1:])+len(head_PYR[1:])+1
resP['AVG'] =resP['SUM']/ntargets
resHe['AVG']=resHe['SUM']/ntargets
resBe['AVG']=resBe['SUM']/ntargets
resC['AVG'] =resC['SUM']/ntargets
resO['AVG'] =resO['SUM']/ntargets

resP['ENERGY']=convMeV(resP['ENERGY'])
resHe['ENERGY']=convMeV(resHe['ENERGY'])/2**(2-alpha)
resBe['ENERGY']=convMeV(resBe['ENERGY'])/4**(2-alpha)
resC['ENERGY']=convMeV(resC['ENERGY'])/6**(2-alpha)
resO['ENERGY']=convMeV(resO['ENERGY'])/8**(2-alpha)

xinterp=resP['ENERGY'].tolist()

resHeint=resHe.copy()
resHeint['ENERGY']=resP['ENERGY']
resHeint['SUM']=interp_XS_EP(xinterp,resHe['ENERGY'],resHe['SUM'],'cubic')
resHeint['AVG']=interp_XS_EP(xinterp,resHe['ENERGY'],resHe['AVG'],'cubic')

resBeint=resBe.copy()
resBeint['ENERGY']=resP['ENERGY']
resBeint['SUM']=interp_XS_EP(xinterp,resBe['ENERGY'],resBe['SUM'],'cubic')
resBeint['AVG']=interp_XS_EP(xinterp,resBe['ENERGY'],resBe['AVG'],'cubic')

resCint=resC.copy()
resCint['ENERGY']=resP['ENERGY']
resCint['SUM']=interp_XS_EP(xinterp,resC['ENERGY'],resC['SUM'],'cubic')
resCint['AVG']=interp_XS_EP(xinterp,resC['ENERGY'],resC['AVG'],'cubic')

resOint=resO.copy()
resOint['ENERGY']=resP['ENERGY']
resOint['SUM']=interp_XS_EP(xinterp,resO['ENERGY'],resO['SUM'],'cubic')
resOint['AVG']=interp_XS_EP(xinterp,resO['ENERGY'],resO['AVG'],'cubic')

EmaxHe=resHe.iloc[-1:]['ENERGY'].tolist()
EmaxBe=resBe.iloc[-1:]['ENERGY'].tolist()
EmaxC=resC.iloc[-1:]['ENERGY'].tolist()
EmaxO=resO.iloc[-1:]['ENERGY'].tolist()

iEmaxHe=resHeint.index[abs(resHeint.loc[:]['ENERGY']-EmaxHe[0])<=0.005].tolist()[-1]
iEmaxBe=resBeint.index[abs(resBeint.loc[:]['ENERGY']-EmaxBe[0])<=0.005].tolist()[-1]
iEmaxC=resCint.index[abs(resCint.loc[:]['ENERGY']-EmaxC[0])<=0.005].tolist()[-1]
iEmaxO=resOint.index[abs(resOint.loc[:]['ENERGY']-EmaxO[0])<=0.005].tolist()[-1]


# plt.plot(resP['ENERGY'][iPe100:],resP['AVG'][iPe100:])
# plt.plot(resHeint['ENERGY'][iPe100:iEmaxHe],resHeint['AVG'][iPe100:iEmaxHe])
# plt.plot(resBeint['ENERGY'][iPe100:iEmaxBe],resBeint['AVG'][iPe100:iEmaxBe])
# plt.plot(resCint['ENERGY'][iPe100:iEmaxC],resCint['AVG'][iPe100:iEmaxC])
# plt.plot(resOint['ENERGY'][iPe100:iEmaxO],resOint['AVG'][iPe100:iEmaxO])
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

xtarg=resP.copy()
xtarg['SUM']=resP['SUM']+resHeint['SUM'][:iEmaxHe]+resBeint['SUM'][:iEmaxBe]+resCint['SUM'][:iEmaxC]+resOint['SUM'][:iEmaxO]
xtarg['AVG']=xtarg['SUM']/(5*ntargets)


ymin=0.02
ymax=2
xmin=0.02
xmax=15
xtext=5
ytext=10
icharge=65/1000
iwidth=0.5
ims=18
imew=0.9

fig = plt.figure(figsize=(16,14))

ax1 = plt.subplot(111)
plt.text(5., 1.2, r'$\alpha=1.2$',fontsize=imolsize+4)
plt.plot(convMeV(PS1int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS1int['SLURACILO'][iPe100:]*convcm2/Zch[0]**alpha/ne_adn[0],
         co_p,linewidth=iwidth, label=r'H$^+$')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS1int['SLURACILO'][iHe100:]*convcm2/Zch[1]**alpha/ne_adn[0],
         co_he,linewidth=iwidth, label=r'He$^{+2}$')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS1int['SLURACILO'][iBe100:]*convcm2/Zch[2]**alpha/ne_adn[0],
         co_be,linewidth=iwidth, label=r'Be$^{+4}$')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS1int['SLURACILO'][iCe100:]*convcm2/Zch[3]**alpha/ne_adn[0],
         co_c,linewidth=iwidth, label=r'C$^{+6}$')
plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS1int['SLURACILO'][iOe100:]*convcm2/Zch[4]**alpha/ne_adn[0],
         co_o,linewidth=iwidth, label=r'O$^{+8}$')

ax1.legend(loc='upper left',bbox_to_anchor=(0.09, 1.09, 0, 0),fontsize=ilegsize,ncol=5)

it=0
for i in head_adn[1:]:
    plt.plot(convMeV(PS1int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS1int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_adn[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS1int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_adn[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS1int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_adn[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS1int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_adn[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS1int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_adn[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

it=0
for i in head_CHn[1:]:
    plt.plot(convMeV(PS2int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS2int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_CHn[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS2int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS2int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_CHn[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS2int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS2int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_CHn[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS2int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS2int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_CHn[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS2int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS2int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_CHn[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

it=0
for i in head_PYR[1:]:
    plt.plot(convMeV(PS3int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS3int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_PYR[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS3int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS3int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_PYR[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS3int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS3int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_PYR[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS3int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS3int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_PYR[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS3int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS3int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_PYR[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

# water
plt.plot(convMeV(PS4int['Energy'][iPe100:])/1**(2-alpha),PS4int['H2O'][iPe100:]*convcm2/1**(alpha)/6,co_p,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS4int['Energy'][iHe100:])/2**(2-alpha),HeS4int['H2O'][iHe100:]*convcm2/2**(alpha)/6,co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS4int['Energy'][iBe100:])/4**(2-alpha),BeS4int['H2O'][iBe100:]*convcm2/4**(alpha)/6,co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS4int['Energy'][iCe100:])/6**(2-alpha),CS4int['H2O'][iCe100:]*convcm2/6**(alpha)/6,co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(OS4int['Energy'][iOe100:])/8**(2-alpha),OS4int['H2O'][iOe100:]*convcm2/8**(alpha)/6,co_o,linewidth=iwidth, label='_nolegend_')

errcte=0.3

# # valores medios por ion
xmean=resP['ENERGY'][iPe100:]
ymean=resP['AVG'][iPe100:]
dymean=ymean*errcte
plt.plot(xmean,ymean,color=co_p)
plt.fill_between(xmean,ymean-dymean,ymean+dymean,color='lightgray')
xmean=resHe['ENERGY'][iHe100:]
ymean=resHe['AVG'][iHe100:]
dymean=ymean*errcte
plt.plot(xmean,ymean,color=co_he)
plt.fill_between(xmean,ymean-dymean,ymean+dymean,color='lightgray')
xmean=resBe['ENERGY'][iBe100:]
ymean=resBe['AVG'][iBe100:]
dymean=ymean*errcte
plt.plot(xmean,ymean,color=co_be)
plt.fill_between(xmean,ymean-dymean,ymean+dymean,color='lightgray')
xmean=resC['ENERGY'][iCe100:]
ymean=resC['AVG'][iCe100:]
dymean=ymean*errcte
plt.plot(xmean,ymean,color=co_c)
plt.fill_between(xmean,ymean-dymean,ymean+dymean,color='lightgray')
xmean=resO['ENERGY'][iOe100:]
ymean=resO['AVG'][iOe100:]
dymean=ymean*errcte
plt.plot(xmean,ymean,color=co_o)
plt.fill_between(xmean,ymean-dymean,ymean+dymean,color='lightgray')

# experiments:

# electrons
Ae, = plt.plot(A_e['Ep'][0:iAe]/1000.,A_e['CS'][0:iAe]/45,color='dimgray',
               marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=6,
               linestyle=' ',label='_nolegend_')
Cite,=plt.plot(C_e['Ep'][iCe:]/1000,C_e['CS'][iCe:]/37,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ',label='e')
Tyme, = plt.plot(T_e['Ep'][0:iTe]/1000,T_e['CS'][0:iTe]/42,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ')
Guane,=plt.plot(G_e['Ep'][0:iGe]/1000,G_e['CS'][0:iGe]/49,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ')
PYRe, = plt.plot(bug17['Ep(keV)']/1000,bug17['e+PY']/28,color='dimgray',
                 marker='>',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
THFe09, = plt.plot(THF_e09['Ep(kev)']/1000,THF_e09['CS']/28,color='dimgray',
                   marker='*',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')
THFe17, = plt.plot(bug17['Ep(keV)']/1000,bug17['e+THF']/28,color='dimgray',
                   marker='<',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')
THFe19, = plt.plot(convMeV(THF_e19['v(au)'][0:iTHFe19]),THF_e19['au'][0:iTHFe19]/28*convcm2,color='dimgray',
                   marker='>',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')

# # gases + H1+
plt.plot(gasesH1['E'],gasesH1['N2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
N2p, = plt.plot(gasesH1['E'],gasesH1['N2'],color=co_p,
                marker=r'$\boxast$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['O2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
O2p, = plt.plot(gasesH1['E'],gasesH1['O2'],color=co_p,
                marker=r'$\boxbox$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['CO'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
COp, = plt.plot(gasesH1['E'],gasesH1['CO'],color=co_p,
                marker=r'$\boxbar$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['CO2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CO2p, = plt.plot(gasesH1['E'],gasesH1['CO2'],color=co_p,
                 marker=r'$\boxbslash$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(CH4_p['E'],CH4_p['TCS']/8,color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4p, = plt.plot(CH4_p['E'],CH4_p['TCS']/8,color=co_p,
                   marker=r'$\boxdot$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# gases + He2+
plt.plot(gasesHe2['E'],gasesHe2['N2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
N2he, = plt.plot(gasesHe2['E'],gasesHe2['N2'],color=co_he,
                 marker=r'$\boxast$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['CO'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
COhe, = plt.plot(gasesHe2['E'],gasesHe2['CO'],color=co_he,
                 marker=r'$\boxbar$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['O2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
O2he, = plt.plot(gasesHe2['E'],gasesHe2['O2'],color=co_he,
                 marker=r'$\boxbox$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['CO2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CO2he, = plt.plot(gasesHe2['E'],gasesHe2['CO2'],color=co_he,
                  marker=r'$\boxbslash$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(CH4_he85['E']/2**(2-alpha),CH4_he85['TCS']/8/2**(alpha),color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4he85, = plt.plot(CH4_he85['E']/2**(2-alpha),CH4_he85['TCS']/8/2**(alpha),color=co_he,
                    marker=r'$\boxdot$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# # methane
plt.plot(CH4_he19['E'],CH4_he19['TCS']/8,color=co_he,
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4he19, = plt.plot(CH4_he19['E'],CH4_he19['TCS']/8,color=co_p,
                    marker=r'$\varowedge$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# water
H2Op86, = plt.plot(H2O_p86['E'],H2O_p86['TCS']/6,color=co_p,
                   marker='p',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',linestyle=' ',label='__nolegend__')
H2Op07, = plt.plot(H2O_p07['E'],H2O_p07['TCS']/6,color=co_p,
                   marker='H',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',linestyle=' ',label='__nolegend__')
plt.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)/6,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Ohe, = plt.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)/6,color=co_he,
                  marker=r'$\odot$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')
plt.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)/6,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oc, = plt.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)/6,color=co_c,
                 marker=r'$\ovee$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')
plt.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)/6,color='white',
         marker='s',markersize=ims-3,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oo, = plt.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)/6,color=co_o,
                 marker=r'$\varobslash$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')

# DNA and RNA
Ap, = plt.plot(A_p['E']/1000,A_p['CS']/45,color=co_p,
               marker='o',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
Up, = plt.plot(U_p['E']/1000,U_p['CS']/36,color=co_p,marker='^',
               markersize=ims,markeredgewidth=2,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
PYRp, = plt.plot(PYR_p['E']/1000,PYR_p['au']/28*convcm2,color=co_p,
                 marker='v',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
THFp, = plt.plot(THF_p['Ep']/1000,THF_p['au']/28*convcm2,color=co_p,
                 marker='D',markersize=ims-2,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
AC, = plt.plot(A_C[0]/6**(2-alpha),A_C[1]/6**(alpha)/45,color=co_c,
               marker='d',markersize=ims,markeredgewidth=imew*2.5,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')


plt.xlabel(r"Impact Energy/Z$^{2-\alpha}$ (MeV/amu)", fontsize=ilabelsize,labelpad=5)
plt.ylabel(r"$\sigma_e/Z^{\alpha}$ ($10^{-16}$ cm$^2$)", fontsize=ilabelsize,labelpad=15)
ax1.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, left=True, right=True);
ax1.tick_params(direction='in',which='minor',length=6,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax) 
plt.yscale('log')
plt.xscale('log')
ax1.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))
# ax1.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

plt.subplots_adjust(left=0.1, right=0.975)

plt.savefig("../error1.eps")
#plt.show()


ymin=0.02
ymax=2
xmin=0.02
xmax=15
xtext=5
ytext=10
icharge=65/1000
iwidth=0.5
ims=18
imew=0.9

fig = plt.figure(figsize=(16,14))

ax1 = plt.subplot(111)
plt.text(5., 1.2, r'$\alpha=1.2$',fontsize=imolsize+4)
plt.plot(convMeV(PS1int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS1int['SLURACILO'][iPe100:]*convcm2/Zch[0]**alpha/ne_adn[0],
         co_p,linewidth=iwidth, label=r'H$^+$')
plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS1int['SLURACILO'][iHe100:]*convcm2/Zch[1]**alpha/ne_adn[0],
         co_he,linewidth=iwidth, label=r'He$^{+2}$')
plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS1int['SLURACILO'][iBe100:]*convcm2/Zch[2]**alpha/ne_adn[0],
         co_be,linewidth=iwidth, label=r'Be$^{+4}$')
plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS1int['SLURACILO'][iCe100:]*convcm2/Zch[3]**alpha/ne_adn[0],
         co_c,linewidth=iwidth, label=r'C$^{+6}$')
plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS1int['SLURACILO'][iOe100:]*convcm2/Zch[4]**alpha/ne_adn[0],
         co_o,linewidth=iwidth, label=r'O$^{+8}$')

ax1.legend(loc='upper left',bbox_to_anchor=(0.09, 1.09, 0, 0),fontsize=ilegsize,ncol=5)

it=0
for i in head_adn[1:]:
    plt.plot(convMeV(PS1int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS1int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_adn[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS1int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS1int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_adn[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS1int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS1int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_adn[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS1int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS1int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_adn[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS1int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS1int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_adn[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

it=0
for i in head_CHn[1:]:
    plt.plot(convMeV(PS2int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS2int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_CHn[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS2int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS2int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_CHn[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS2int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS2int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_CHn[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS2int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS2int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_CHn[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS2int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS2int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_CHn[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

it=0
for i in head_PYR[1:]:
    plt.plot(convMeV(PS3int['ENERGY'][iPe100:])/Zch[0]**(2-alpha),PS3int[i][iPe100:]*convcm2/Zch[0]**alpha/ne_PYR[it],
             co_p,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(HeS3int['ENERGY'][iHe100:])/Zch[1]**(2-alpha),HeS3int[i][iHe100:]*convcm2/Zch[1]**alpha/ne_PYR[it],
             co_he,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(BeS3int['ENERGY'][iBe100:])/Zch[2]**(2-alpha),BeS3int[i][iBe100:]*convcm2/Zch[2]**alpha/ne_PYR[it],
             co_be,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(CS3int['ENERGY'][iCe100:])/Zch[3]**(2-alpha),CS3int[i][iCe100:]*convcm2/Zch[3]**alpha/ne_PYR[it],
             co_c,linewidth=iwidth, label='_nolegend_')
    plt.plot(convMeV(OS3int['ENERGY'][iOe100:])/Zch[4]**(2-alpha),OS3int[i][iOe100:]*convcm2/Zch[4]**alpha/ne_PYR[it],
             co_o,linewidth=iwidth, label='_nolegend_')
    it=it+1

# water
plt.plot(convMeV(PS4int['Energy'][iPe100:])/1**(2-alpha),PS4int['H2O'][iPe100:]*convcm2/1**(alpha)/6,co_p,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(HeS4int['Energy'][iHe100:])/2**(2-alpha),HeS4int['H2O'][iHe100:]*convcm2/2**(alpha)/6,co_he,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(BeS4int['Energy'][iBe100:])/4**(2-alpha),BeS4int['H2O'][iBe100:]*convcm2/4**(alpha)/6,co_be,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(CS4int['Energy'][iCe100:])/6**(2-alpha),CS4int['H2O'][iCe100:]*convcm2/6**(alpha)/6,co_c,linewidth=iwidth, label='_nolegend_')
plt.plot(convMeV(OS4int['Energy'][iOe100:])/8**(2-alpha),OS4int['H2O'][iOe100:]*convcm2/8**(alpha)/6,co_o,linewidth=iwidth, label='_nolegend_')

# # valores medios por ion
# plt.plot(resP['ENERGY'][iPe100:],resP['AVG'][iPe100:],color=co_p)
# plt.plot(resHe['ENERGY'][iHe100:],resHe['AVG'][iHe100:],color=co_he)
# plt.plot(resBe['ENERGY'][iBe100:],resBe['AVG'][iBe100:],color=co_be)
# plt.plot(resC['ENERGY'][iCe100:],resC['AVG'][iCe100:],color=co_c)
# plt.plot(resO['ENERGY'][iOe100:],resO['AVG'][iOe100:],color=co_o)

errcte=0.4

# valor medio total
xmean=resP['ENERGY'][iPe100:]
ymean=xtarg['AVG'][iPe100:]
dymean=xtarg['AVG'][iPe100:]*errcte
plt.plot(xmean,ymean,color='k')
plt.fill_between(xmean,ymean-dymean,ymean+dymean,color='lightgray')

# experiments:

# electrons
Ae, = plt.plot(A_e['Ep'][0:iAe]/1000.,A_e['CS'][0:iAe]/45,color='dimgray',
               marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=6,
               linestyle=' ',label='_nolegend_')
Cite,=plt.plot(C_e['Ep'][iCe:]/1000,C_e['CS'][iCe:]/37,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ',label='e')
Tyme, = plt.plot(T_e['Ep'][0:iTe]/1000,T_e['CS'][0:iTe]/42,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ')
Guane,=plt.plot(G_e['Ep'][0:iGe]/1000,G_e['CS'][0:iGe]/49,color='dimgray',
         marker='s',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',markevery=4,
         linestyle=' ')
PYRe, = plt.plot(bug17['Ep(keV)']/1000,bug17['e+PY']/28,color='dimgray',
                 marker='>',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
THFe09, = plt.plot(THF_e09['Ep(kev)']/1000,THF_e09['CS']/28,color='dimgray',
                   marker='*',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')
THFe17, = plt.plot(bug17['Ep(keV)']/1000,bug17['e+THF']/28,color='dimgray',
                   marker='<',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')
THFe19, = plt.plot(convMeV(THF_e19['v(au)'][0:iTHFe19]),THF_e19['au'][0:iTHFe19]/28*convcm2,color='dimgray',
                   marker='>',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                   linestyle=' ',label='_nolegend_')

# gases + H1+
plt.plot(gasesH1['E'],gasesH1['N2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
N2p, = plt.plot(gasesH1['E'],gasesH1['N2'],color=co_p,
                marker=r'$\boxast$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['O2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
O2p, = plt.plot(gasesH1['E'],gasesH1['O2'],color=co_p,
                marker=r'$\boxbox$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['CO'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
COp, = plt.plot(gasesH1['E'],gasesH1['CO'],color=co_p,
                marker=r'$\boxbar$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesH1['E'],gasesH1['CO2'],color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CO2p, = plt.plot(gasesH1['E'],gasesH1['CO2'],color=co_p,
                 marker=r'$\boxbslash$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(CH4_p['E'],CH4_p['TCS']/8,color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4p, = plt.plot(CH4_p['E'],CH4_p['TCS']/8,color=co_p,
                   marker=r'$\boxdot$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# gases + He2+
plt.plot(gasesHe2['E'],gasesHe2['N2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
N2he, = plt.plot(gasesHe2['E'],gasesHe2['N2'],color=co_he,
                 marker=r'$\boxast$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['CO'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
COhe, = plt.plot(gasesHe2['E'],gasesHe2['CO'],color=co_he,
                 marker=r'$\boxbar$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['O2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
O2he, = plt.plot(gasesHe2['E'],gasesHe2['O2'],color=co_he,
                 marker=r'$\boxbox$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(gasesHe2['E'],gasesHe2['CO2'],color='white',
         marker='s',markersize=ims-1,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CO2he, = plt.plot(gasesHe2['E'],gasesHe2['CO2'],color=co_he,
                  marker=r'$\boxbslash$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
plt.plot(CH4_he85['E']/2**(2-alpha),CH4_he85['TCS']/8/2**(alpha),color='white',
        marker='s',markersize=ims-2,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4he85, = plt.plot(CH4_he85['E']/2**(2-alpha),CH4_he85['TCS']/8/2**(alpha),color=co_he,
                    marker=r'$\boxdot$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# # methane
plt.plot(CH4_he19['E'],CH4_he19['TCS']/8,color=co_he,
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')
CH4he19, = plt.plot(CH4_he19['E'],CH4_he19['TCS']/8,color=co_p,
                    marker=r'$\varowedge$',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ')

# water
H2Op86, = plt.plot(H2O_p86['E'],H2O_p86['TCS']/6,color=co_p,
                   marker='p',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',linestyle=' ',label='__nolegend__')
H2Op07, = plt.plot(H2O_p07['E'],H2O_p07['TCS']/6,color=co_p,
                   marker='H',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',linestyle=' ',label='__nolegend__')
plt.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)/6,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Ohe, = plt.plot(H2O_he2['E']/2**(2-alpha),H2O_he2['H2O']/2**(alpha)/6,color=co_he,
                  marker=r'$\odot$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')
plt.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)/6,color='white',
         marker='o',markersize=ims,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oc, = plt.plot(H2O_c6['E']/6**(2-alpha),H2O_c6['TCS']/6**(alpha)/6,color=co_c,
                 marker=r'$\ovee$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')
plt.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)/6,color='white',
         marker='s',markersize=ims-3,markeredgewidth=imew,markerfacecolor='white',linestyle=' ',label='__none__')
H2Oo, = plt.plot(H2O_o8['E']/8**(2-alpha),H2O_o8['TCS']/8**(alpha)/6,color=co_o,
                 marker=r'$\varobslash$',markersize=ims,markeredgewidth=imew*2.0,linestyle=' ',label='__nolegend__')

# DNA and RNA
Ap, = plt.plot(A_p['E']/1000,A_p['CS']/45,color=co_p,
               marker='o',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
Up, = plt.plot(U_p['E']/1000,U_p['CS']/36,color=co_p,marker='^',
               markersize=ims,markeredgewidth=2,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')
PYRp, = plt.plot(PYR_p['E']/1000,PYR_p['au']/28*convcm2,color=co_p,
                 marker='v',markersize=ims,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
THFp, = plt.plot(THF_p['Ep']/1000,THF_p['au']/28*convcm2,color=co_p,
                 marker='D',markersize=ims-2,markeredgewidth=imew*2.0,markerfacecolor='white',
                 linestyle=' ',label='_nolegend_')
AC, = plt.plot(A_C[0]/6**(2-alpha),A_C[1]/6**(alpha)/45,color=co_c,
               marker='d',markersize=ims,markeredgewidth=imew*2.5,markerfacecolor='white',
               linestyle=' ',label='_nolegend_')


plt.xlabel(r"Impact Energy/Z$^{2-\alpha}$ (MeV/amu)", fontsize=ilabelsize,labelpad=5)
plt.ylabel(r"$\sigma_e/Z^{\alpha}$ ($10^{-16}$ cm$^2$)", fontsize=ilabelsize,labelpad=15)
ax1.tick_params(direction='in',which='major',labelsize=itickssize,length=8,bottom=True, top=True, left=True, right=True);
ax1.tick_params(direction='in',which='minor',length=6,bottom=True, top=True, left=True, right=True);
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax) 
plt.yscale('log')
plt.xscale('log')
ax1.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))
# ax1.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,g}'))

plt.subplots_adjust(left=0.1, right=0.975)

plt.savefig("../error2.eps")
#plt.show()

