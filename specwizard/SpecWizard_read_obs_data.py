
import h5py 
import numpy as np
import matplotlib.pyplot as plt


class read_obs_data():
    '''
    This class reads the data from observational data from Kim et al 2020 and Vid Iršič et al 2017
    '''
    def __init__(self, datapath='./data/'):
        #
        self.datapath = datapath
        
    def kim2020(self,filename='kim2020_data',plot=False):
        kim2020_hdf5 = h5py.File(self.datapath+filename+".hdf5",'r')
        kim2020_dic = {}
        
        #reads table1 
        kim2020_dic['Table1']         = {}
        kim2020_dic['Table1']['Info'] = kim2020_hdf5['Table1'].attrs['Info']
        kim2020_dic['Table1']['AGN']  = {}
        list_agn                      = kim2020_hdf5['Table1']['AGN'].keys()
        for agn in list_agn:
            kim2020_dic['Table1']['AGN'][agn] = {}
            list_attributes = list(kim2020_hdf5['Table1']['AGN'][agn].attrs)
            for att in list_attributes:
                    kim2020_dic['Table1']['AGN'][agn][att] = kim2020_hdf5['Table1']['AGN'][agn].attrs[att]
                    
        #reads table4 
        kim2020_dic['Table4'] = {}
        kim2020_dic['Table4']['Info'] = kim2020_hdf5['Table4'].attrs['Info']
        kim2020_dic['Table4']['Redshifts'] = {}
        list_z   = kim2020_hdf5['Table4']['Redshifts'].keys()
        for zz in list_z:
            kim2020_dic['Table4']['Redshifts'][zz] = {}
            list_attributes = list(kim2020_hdf5['Table4']['Redshifts'][zz].attrs)
            for att in list_attributes:
                    kim2020_dic['Table4']['Redshifts'][zz][att] = kim2020_hdf5['Table4']['Redshifts'][zz].attrs[att]
        
        #reads table5
        kim2020_dic['Table5'] = {}
        kim2020_dic['Table5']['Info'] = kim2020_hdf5['Table5'].attrs['Info']
        kim2020_dic['Table5']['Flux'] = {}
        kim2020_dic['Table5']['Flux']['Value'] = kim2020_hdf5['Table5']['Flux'][...]
        kim2020_dic['Table5']['Flux']['Info'] = kim2020_hdf5['Table5']['Flux'].attrs['Info']

        kim2020_dic['Table5']['Redshifts'] = {}
        list_z = kim2020_hdf5['Table5']['Redshifts'].keys()
        for zz in list_z:
            kim2020_dic['Table5']['Redshifts'][zz] = {}
            list_attributes = list(kim2020_hdf5['Table5']['Redshifts'][zz].attrs)
            for att in list_attributes:
                    kim2020_dic['Table5']['Redshifts'][zz][att] = kim2020_hdf5['Table5']['Redshifts'][zz].attrs[att]
            for ks in kim2020_hdf5['Table5']['Redshifts'][zz].keys():
                kim2020_dic['Table5']['Redshifts'][zz][ks] = {}
                kim2020_dic['Table5']['Redshifts'][zz][ks]['Value'] =  kim2020_hdf5['Table5']['Redshifts'][zz][ks][...]
                kim2020_dic['Table5']['Redshifts'][zz][ks]['Info']  =  kim2020_hdf5['Table5']['Redshifts'][zz][ks].attrs['Info']
        if plot==True:
            fig, axs = plt.subplots(1,7, figsize=(20, 3), facecolor='w', edgecolor='k',sharey=True)
            fig.subplots_adjust(hspace = 1.5, wspace=0)
            fig.text(0.095, 0.5, '<P(F)>', va='center', rotation='vertical')

            redshifts = kim2020_dic['Table5']['Redshifts'].keys()
            flux_bins = kim2020_dic['Table5']['Flux']['Value']

            for ax,zz in zip(axs,redshifts):
                P_f = kim2020_dic['Table5']['Redshifts'][zz]['P(F)']['Value']
                pferror = kim2020_dic['Table5']['Redshifts'][zz]['error']['Value']
                ax.errorbar(flux_bins,P_f, yerr=pferror, fmt="o",c='r',ms=2.5,label='Kim 2012')

                ax.set_title("z: "+zz)
                ax.set_yscale('log')
                ax.set_ylim(1e-3,6)
                ax.set_xlim(0,1)
                ax.set_xlabel("F")
                x_ticks = ax.xaxis.get_major_ticks()
                x_ticks[-1].label1.set_visible(False)
                if zz == '3.36':
                    ax.legend()
        
        return kim2020_dic 
    def Day19(self,filename='Day19.csv'):

        # power spectrum at z=3.13 from Day19

        filename  = self.datapath+filename
        data      = np.loadtxt(filename, delimiter=",", comments='#', ).T
        x         = np.array(data[0,:])
        y         = np.array(data[1,:])
        nx        = int(len(x)/2)
        logk      = x[0:nx]
        Pmin      = y[0:nx]
        logkm     = x[nx:]
        Pmax      = y[nx:]
        Pmean     = 0.5*(Pmin+Pmax)
        LogVar    = logk + Pmean - np.log10(np.pi)
        Min       = logk + Pmin  - np.log10(np.pi)
        Max       = logk + Pmax  - np.log10(np.pi)
        info      = 'Power spectrum at z=3.13 from Day19'
        Day19     = {'logk':logk, 'logVar':LogVar, 'Min':Min, 'Max':Max, 'LogP': Pmean, 'LogPmin':Pmin, 'LogPmax':Pmax,'info':info}
        return Day19

    def Becker2014(self,filename='Becker2014.csv'):
        becker2014 = np.genfromtxt(self.datapath+filename, delimiter=',')

        info = 'Average flux data from Becker et al 2014'
        beck_dic = {'z':becker2014[:,0],'F':becker2014[:,1],'info':info}

        return beck_dic
    def Walther18(self,filename='Walther18.csv'):
        ## power spectrum from Walther at z=3.0
        filename  = self.datapath+filename

        data      = np.loadtxt(filename, delimiter=",", comments='#', ).T
        x         = np.array(data[0,:])
        y         = np.array(data[1,:])
        nx        = int(len(x)/2)
        logk      = x[0:nx]
        Pmin      = y[0:nx]
        logkm     = x[nx:]
        Pmax      = y[nx:]
        Pmean     = 0.5*(Pmin+Pmax)
        info      = ' power spectrum from Walther at z=3.0'
        Walther18 = {'logk':logk, 'logVar':Pmean, 'Min':Pmin, 'Max':Pmax,'info':info}
        return Walther18

    def VidIrsic_2017(self,filename='VidIrsic2017_data',plot=False):
        
        VidIrsic2017_data = h5py.File(self.datapath+filename+".hdf5",'r')
        VidIrsic_2017_dic = {}
        
        VidIrsic_2017_dic['Figure6'] = {}
        VidIrsic_2017_dic['Figure6']['Info']     = VidIrsic2017_data['Figure6'].attrs['Info']
        VidIrsic_2017_dic['Figure6']['Redshifts'] = {}
        list_z = VidIrsic2017_data['Figure6']['Redshifts'].keys()
        for zz in list_z:
            VidIrsic_2017_dic['Figure6']['Redshifts'][zz] = {}
            Lisattr = list(VidIrsic2017_data['Figure6']['Redshifts'][zz].attrs)
            for latr in Lisattr:
                VidIrsic_2017_dic['Figure6']['Redshifts'][zz][latr] = VidIrsic2017_data['Figure6']['Redshifts'][zz].attrs[latr]
       
        VidIrsic_2017_dic['TableA1'] = {}
        VidIrsic_2017_dic['TableA1']['Info'] = VidIrsic2017_data['TableA1'].attrs['Info']
        VidIrsic_2017_dic['TableA1']['Redshifts'] = {}
        list_z = VidIrsic2017_data['TableA1']['Redshifts'].keys()
        for zz in list_z:
            VidIrsic_2017_dic['TableA1']['Redshifts'][zz] = {}
            list_keys = VidIrsic2017_data['TableA1']['Redshifts'][zz].keys()
            for datafield in list_keys:
                VidIrsic_2017_dic['TableA1']['Redshifts'][zz][datafield] = {}
                VidIrsic_2017_dic['TableA1']['Redshifts'][zz][datafield]['Value'] = VidIrsic2017_data['TableA1']['Redshifts'][zz][datafield][...]
                list_attributes = list(VidIrsic2017_data['TableA1']['Redshifts'][zz][datafield].attrs)
                for att in list_attributes:
                    VidIrsic_2017_dic['TableA1']['Redshifts'][zz][datafield][att] = VidIrsic2017_data['TableA1']['Redshifts'][zz][datafield].attrs[att]
        
        if plot==True:
            #redshifts = VidIrsic_2017_dic['TableA1']['Redshifts'].keys()
            redshifts = ['3.0','3.6','4.2']
            for zz in redshifts:
                kvals = VidIrsic_2017_dic['TableA1']['Redshifts'][zz]['k']['Value']
                PFvals = VidIrsic_2017_dic['TableA1']['Redshifts'][zz]['P_F']['Value']
                errstat = VidIrsic_2017_dic['TableA1']['Redshifts'][zz]['sigma_stat']['Value']
                errsys = VidIrsic_2017_dic['TableA1']['Redshifts'][zz]['sigma_sys']['Value']
                errtot = np.sqrt(errstat**2+errsys**2)

                plt.errorbar(kvals,(kvals*PFvals)/np.pi,(kvals*errtot)/np.pi,fmt = '-o',ms=4,label='z: '+zz)
            plt.legend()
            plt.ylabel(r'$kP_F(k)/\pi$')
            plt.xlabel(r'$k$')
            plt.ylim(0.01,0.3)
        return VidIrsic_2017_dic