#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Tkinter as tk
import ttk
import matplotlib.pyplot as plt
import h5py
import os
import scipy.signal as signal
from scipy import sparse
from scipy.sparse.linalg import spsolve
import tkFileDialog
import glob
import tkMessageBox
import numpy as np
import shutil
from decimal import Decimal
from PIL import Image


# In[6]:


LARGE_FONT=('Verdana', 16)
color='blue'
small_font=('Verdana', 10)
number_font=('Verdana', 6)


# In[7]:


def elliptic_bandpass(order, rp, rs, lowcut, highcut):
    fs=256
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.ellip(order, rp, rs, [low, high], btype='bandpass', analog=False, output='ba')
    return b, a



def power(x):
    x=1.0*(x)
    return sum(x**2)/len(x)

def CAR(samples):
    s=[]
    for i in range(0, len(samples)):
        s.append((sum(samples))/64.0)
    
    sample=samples-s
    return sample, s




def baseline_als(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    #D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z


spr_channels=[0,2,7,9,11,13,15,62,25,27,29,31,33,63,43,45,47,49,51,57,59]
spr_channels=np.array(spr_channels)
spr_channel_names=['fp1','fp2','f7','f3','fz','f4','f8','A1','T3','C3','Cz','C4','T4','A2','T5','P3','Pz','P4','T6'
                  ,'O1','O2']


# In[18]:


class Qeeg(tk.Tk):

    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "QUANTITATIVE EEG")


        root = tk.Frame(self, relief='raised', borderwidth=5)
        root.pack(side="top", fill="both", expand = True)
        root.grid_rowconfigure(1, weight=1)
        root.grid_columnconfigure(1, weight=1)
        #self.dirs={1:'N.A.', 2:'N.A.', 7:'N.A.', 30:'N.A.', 90:'N.A.'}
        self.frames={}
        for F in (StartPage, QeegPage):

            frame=F(root, self)
            self.frames[F]= frame
            frame.grid(row=0, column=0, sticky='nsew')

        self.show_frame(StartPage)
        self.filepath = tk.StringVar()
        

    def show_frame(self, cont):
        #self.dirs={1:'N.A.', 2:'N.A.', 7:'N.A.', 30:'N.A.', 90:'N.A.'}
        frame=self.frames[cont]
        frame.tkraise()
        

    def quit(self, frame):
        app.destroy()

class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        
        label = tk.Label(self, text="Choose the Folder", font=LARGE_FONT)
        label.grid(row=1, column=50, sticky='N')
        self.dirs={1:'N.A.', 2:'N.A.', 7:'N.A.', 30:'N.A.', 90:'N.A.'}
        button = ttk.Button(self, text="Browse Folder",
                            command=lambda: self.browse_folder(controller, 0))

        button.grid(row=2, column=50)
        button = ttk.Button(self, text="Browse",
                            command=lambda: self.browse_folder(controller, 1))

        button.grid(row=3, column=51)
        button = ttk.Button(self, text="Browse",
                            command=lambda: self.browse_folder(controller, 2))

        button.grid(row=4, column=51)
        button = ttk.Button(self, text="Browse",
                            command=lambda: self.browse_folder(controller, 7))

        button.grid(row=5, column=51)
        button = ttk.Button(self, text="Browse",
                            command=lambda: self.browse_folder(controller, 30))

        button.grid(row=6, column=51)
        button = ttk.Button(self, text="Browse",
                            command=lambda: self.browse_folder(controller, 90))

        button.grid(row=7, column=51)
        
        
        
        #label= tk.Label(self, text="Select the Folder with patient's data for all days")
        #label.grid(row=0, column=0, sticky='NW')
        #self.grid_rowconfigure(1, weight=1)
        #self.grid_columnconfigure(1, weight=1)
        c=2
        for i in [1, 2, 7, 30, 90]:
            c=c+1
            label= tk.Label(self, text='DAY'+ str(i)+':')
            label.grid(row= c, column=30, sticky='NSEW')

        self.filepath = tk.StringVar()
        self.e1= tk.Entry(self, text='Day 1 Folder')
        self.e2= tk.Entry(self, text='Day 2 Folder')
        self.e3= tk.Entry(self, text='Day 7 Folder')
        self.e4= tk.Entry(self, text='Day 30 Folder')
        self.e5= tk.Entry(self, text='Day 90 Folder')

        self.e1.grid(row= 3, column=50, sticky='NSEW')
        self.e2.grid(row= 4, column=50, sticky='NSEW')
        self.e3.grid(row= 5, column=50, sticky='NSEW')
        self.e4.grid(row= 6, column=50, sticky='NSEW')
        self.e5.grid(row= 7, column=50,sticky='NSEW')

        #self.dirs={}

        button2 = ttk.Button(self, text="Compute QEEG", command=lambda: self.compute(controller))
        button2.grid(row= 9, column= 50)
        
        #self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)
        self.grid_rowconfigure(4, weight=1)
        self.grid_rowconfigure(5, weight=1)
        self.grid_rowconfigure(6, weight=1)
        self.grid_rowconfigure(7, weight=1)
        self.grid_rowconfigure(8, weight=1)
        self.grid_columnconfigure(50, weight=1)
        #self.grid_columnconfigure(1, weight=1)

    def browse_folder(self, controller, br):
        
        if br==0:
            self.filepath=tkFileDialog.askdirectory()
            arr= glob.glob(self.filepath+'/*/')
            arr.sort()
            self.dirs={1:'N.A.', 2:'N.A.', 7:'N.A.', 30:'N.A.', 90:'N.A.'}
            if arr!=None:
                for f in arr:
                    f1=f.lower()
                    if f1.find('day1') != -1 :

                        self.dirs[1]=f

                    elif f1.find('day2') != -1:

                        self.dirs[2]=f

                    elif f1.find('day7') != -1 :

                        self.dirs[7]=f

                    elif f1.find('day30') != -1 :

                        self.dirs[30]=f

                    elif f1.find('day90') != -1 :

                        self.dirs[90]=f
                        
                

            if (self.dirs[1]=='N.A.' and self.dirs[2]=='N.A.' and self.dirs[7]=='N.A.' and 
                self.dirs[30]=='N.A.' and self.dirs[90]=='N.A.'):
                tkMessageBox.showerror("Error", 'NO Directory of Patient Data Found')
                
        else:
            self.filepath=tkFileDialog.askdirectory()
            self.dirs[br]=self.filepath
            if br==1:
                self.e1.delete(0, 'end')
            elif br==2:
                self.e1.delete(0, 'end')
            elif br==3:
                self.e1.delete(0, 'end')
                
            elif br==4:
                self.e1.delete(0, 'end')
            else:
                self.e5.delete(0, 'end')
        
        controller.dirs=self.dirs
        controller.filepath=self.filepath
        self.e1.delete(0, 'end')
        self.e2.delete(0, 'end')
        self.e3.delete(0, 'end')
        self.e4.delete(0, 'end')
        self.e5.delete(0, 'end')
        self.e1.insert(0, self.dirs[1])
        self.e2.insert(0, self.dirs[2])
        self.e3.insert(0, self.dirs[7])
        self.e4.insert(0, self.dirs[30])
        self.e5.insert(0, self.dirs[90])
        
        print controller.dirs
        controller.dirs=self.dirs
        controller.filepath=self.filepath
        
       
    
        
    def elliptic_bandpass(self, order, rp, rs, lowcut, highcut):
        fs=256
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = signal.ellip(order, rp, rs, [low, high], btype='bandpass', analog=False, output='ba')
        return b, a



    def power(self, x):
        x=1.0*(x)
        return sum(x**2)/len(x)
        
    def subpower(self, efile, spr_channel):
        spr_channels=[0,2,7,9,11,13,15,62,25,27,29,31,33,63,43,45,47,49,51,57,59]
        spr_channels=np.array(spr_channels)
        powerr=[]
        for i in range (0, len(efile[0]), 2):
            activity=str(efile[0][i][72:])
            eeg=[]
            for j in range(0, len(efile)):
                eeg.append(h5py.File(efile[j][i],mode='r'))
                eeg.append(h5py.File(efile[j][i+1],mode='r'))
            samples=[]
            sampletime1=[]
            
            for j in range(0, len(eeg)):
                sam=[]
                sa=np.array(eeg[j]["RawData"]['Samples'])
                sa=sa.T
                #sa,_=CAR(sa)

                #for i in range(0, len(sa)):
                 #   z=baseline_als(sa[i], 10**4, 0.1, 5)
                  #  sa[i]=sa[i]-z

                samples.append(sa)
                sampletime = np.array(eeg[j]['AsynchronData']["Time"])
                sam.append(int(sampletime[0]))
                sam.append(int(sampletime[1]))
                sampletime1.append(sam)

            freq=[[0.5, 4], [4, 8], [8, 13], [13, 30], [30, 100]]
            fs=256
            order=4
            rp=0.5
            rs=30
            powe=[]
            for j in range(0, len(samples)):

                alpha=[]
                beta=[]
                gamma=[]
                delta=[]
                theta=[]
                p=[]
                k = spr_channel
                b, a = self.elliptic_bandpass(order, rp, rs, freq[0][0], freq[0][1])
                delta.append(signal.filtfilt(b, a, samples[j][k], padlen=0))
                b, a = self.elliptic_bandpass(order, rp, rs, freq[1][0], freq[1][1])
                theta.append(signal.filtfilt(b, a, samples[j][k], padlen=0))
                b, a = self.elliptic_bandpass(order, rp, rs, freq[2][0], freq[2][1])
                alpha.append(signal.filtfilt(b, a, samples[j][k], padlen=0))
                b, a = self.elliptic_bandpass(order, rp, rs, freq[3][0], freq[3][1])
                beta.append(signal.filtfilt(b, a, samples[j][k], padlen=0))
                b, a = self.elliptic_bandpass(order, rp, rs, freq[4][0], freq[4][1])
                gamma.append(signal.filtfilt(b, a, samples[j][k], padlen=0))

                k=0
                alpha[k]=alpha[k][sampletime1[j][0]:sampletime1[j][1]]
                beta[k]=beta[k][sampletime1[j][0]:sampletime1[j][1]]
                gamma[k]=gamma[k][sampletime1[j][0]:sampletime1[j][1]]
                delta[k]=delta[k][sampletime1[j][0]:sampletime1[j][1]]
                theta[k]=theta[k][sampletime1[j][0]:sampletime1[j][1]]



                p.append(self.power(delta[k]))
                p.append(self.power(theta[k]))
                p.append(self.power(alpha[k]))
                p.append(self.power(beta[k]))
                p.append(self.power(gamma[k]))

                powe.append(p)
            powerr.append(powe)
        
        
        
        return powerr
    
    def compute_bsi_t(self, bsir, t, T):
        bsi_channels=[[13, 31], [9, 27], [31, 49], [27, 45], [49, 59], [45, 57], [13, 33], [9, 25]]
        ni=np.zeros(T)
        for i in range(0, len(bsi_channels)):

            rj=np.abs(np.fft.fft(bsir[bsi_channels[i][0]][t-T:t]))
            lj=np.abs(np.fft.fft(bsir[bsi_channels[i][1]][t-T:t]))
            ni=ni+((rj-lj)/(rj+lj))/len(bsi_channels)


        ni=np.mean(np.abs(ni))
        #print (ni.real)
        return ni
    
    def compute_bsi(self, controller, efile):
        bsi_channels=[[13, 31], [9, 27], [31, 49], [27, 45], [49, 59], [45, 57], [13, 33], [9, 25]]
        bsi_channel_names=[['F4', 'C4'], ['F3', 'C3'], ['C4', 'P4'], ['C3', 'P3'], ['P4', 'O2'], ['P3', 'O1'], ['F4', 'T4'], 
                   ['F3', 'T3']]
        bsifile=[]
        for i in range(0, len(efile)):
            eeg=h5py.File(efile[i], mode='r')
            samples=np.array(eeg['RawData']['Samples'])
            samples=samples.T
            sampletime = np.array(eeg['AsynchronData']["Time"])

            bsir=[]
            fs=256
            order=4
            rp=0.5
            rs=30
            for i in range (0, len(samples)):
                b, a= self.elliptic_bandpass(order, rp, rs, 1, 25)
                bsir.append(signal.filtfilt(b, a, samples[i], padlen=0))
                
            
            T=256*4  #window time
            bsi=[]
            for i in range(sampletime[0]+T, sampletime[1], 512):
                bsi.append(self.compute_bsi_t(bsir, i, T))
                
            bsifile.append(bsi)
            
        return bsifile

    def compute(self, controller):
        ##code to comptue all_powers
        efile=[]
        n=0
        for i in sorted(self.dirs):
            if self.dirs[i]!= 'N.A.':
                arr=list(glob.glob(self.dirs[i]+'/*hdf5'))
                arr.sort()
                n=n+1
                efile.append(arr[len(arr)-14:])
                
        if n==0:
            tkMessageBox.showerror("Error", 'NO Directory of Patient Data Found')
        spr_channels=[0,2,7,9,11,13,15,62,25,27,29,31,33,63,43,45,47,49,51,57,59]
        spr_channels=np.array(spr_channels)
        all_power=[]
        progressbar=ttk.Progressbar(self, orient='horizontal', length=300, mode='determinate')
        progressbar.grid(row=10, column=50, sticky='WE')
        progressbar['value']=0
        progressbar['maximum']=len(spr_channels)
        divisions=len(spr_channels)
        for i in range(0, len(spr_channels)):
            #print 1
            powerr=self.subpower(efile, i)
            all_power.append(powerr)
            progressbar['value']=i
            progressbar.update()

        bsifile=[]
        for i in range(0, len(efile)):
            bsifile.append(self.compute_bsi(controller, efile[i]))
            progressbar['value']=i
            progressbar.update()
        
        
        
        folder=(self.filepath+'/qeeg')
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)
        np.save(folder+'/all_powers', all_power)
        np.save(folder+'/bsi', bsifile)
        controller.show_frame(QeegPage)

class QeegPage(tk.Frame):

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)
        #label = tk.Label(self, text="Page One!!!", font=LARGE_FONT)
        #label.pack(pady=10,padx=10)
        button0= ttk.Button(self, text="Take Data\nofanotherPatient",
                            command=lambda: self.show_another(controller))
        button0.grid(row=0, column=0, sticky='NS')
        button1 = ttk.Button(self, text="Show Results\nEYE Closed(a)",
                            command=lambda: self.show_results_a(controller, 0))
        button1.grid(row=0, column=18, sticky='NS')
        button2 = ttk.Button(self, text="Show Results\nEYE Open(b)",
                            command=lambda: self.show_results_a(controller, 1))
        button2.grid(row=0, column=19, sticky='NS')
        
        
        self.dayframes=[]
        days=[1, 2, 7, 30, 90]
        for i in range(0, 5):
            dayframe= tk.LabelFrame(self, text=' Day '+str(days[i]), font=LARGE_FONT)
            dayframe.grid(row=1, column=i*9, columnspan=8)
            
            inFileLbl = tk.Label(dayframe, text='Frequency\n(Hz)',fg=color, font=small_font)
            inFileLbl.grid(row=0, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Right',fg=color, font=small_font)
            inFileLbl.grid(row=0, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Hemisphere',fg=color, font=small_font)
            inFileLbl.grid(row=0, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Pre')
            inFileLbl.grid(row=1, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Post')
            inFileLbl.grid(row=1, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Left',fg=color, font=small_font)
            inFileLbl.grid(row=0, column=4, sticky='WE', padx=0, pady=2)            
            inFileLbl = tk.Label(dayframe, text='Hemisphere',fg=color, font=small_font)
            inFileLbl.grid(row=0, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Pre')
            inFileLbl.grid(row=1, column=4, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Post')
            inFileLbl.grid(row=1, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Mean\nDelta')
            inFileLbl.grid(row=3, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Mean\nTheta')
            inFileLbl.grid(row=6, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Mean\nAlpha')
            inFileLbl.grid(row=10, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Mean\nBeta')
            inFileLbl.grid(row=14, column=0, sticky='WE', padx=0, pady=2)
            
            r=20
            inFileLbl = tk.Label(dayframe, text='Time\nCompressed', fg=color, font=small_font)
            inFileLbl.grid(row=r+0, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Right',fg=color, font=small_font)
            inFileLbl.grid(row=r+0, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Hemisphere',fg=color, font=small_font)
            inFileLbl.grid(row=r+0, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Pre')
            inFileLbl.grid(row=r+1, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Post')
            inFileLbl.grid(row=r+1, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Left',fg=color, font=small_font)
            inFileLbl.grid(row=r+0, column=4, sticky='WE', padx=0, pady=2)            
            inFileLbl = tk.Label(dayframe, text='Hemisphere',fg=color, font=small_font)
            inFileLbl.grid(row=r+0, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Pre')
            inFileLbl.grid(row=r+1, column=4, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Post')
            inFileLbl.grid(row=r+1, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Delta\n%')
            inFileLbl.grid(row=r+3, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Theta\n%')
            inFileLbl.grid(row=r+6, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Alpha\n%')
            inFileLbl.grid(row=r+10, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Beta\n%')
            inFileLbl.grid(row=r+14, column=0, sticky='WE', padx=0, pady=2)
            
            r1=40
            inFileLbl = tk.Label(dayframe, text='Indices',fg=color, font= small_font)
            inFileLbl.grid(row=r1+0, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Power Ratio\nIndex')
            inFileLbl.grid(row=r1+3, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Delta/Alpha\nRatio')
            inFileLbl.grid(row=r1+6, column=0, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Pre')
            inFileLbl.grid(row=r1+8, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Post')
            inFileLbl.grid(row=r1+8, column=4, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='Brain\nSymmetry\nIndex')
            inFileLbl.grid(row=r1+10, column=0, sticky='WE', padx=0, pady=2)
            
            self.dayframes.append(dayframe)
       # print controller.dirs
    
    def show_another(self, controller):
        #controller.dirs={1:'N.A.', 2:'N.A.', 7:'N.A.', 30:'N.A.', 90:'N.A.'}
        controller.show_frame(StartPage)

    def show_results_a(self, controller, ab):
        key=controller.dirs.keys()
        key.sort()
        keys=[]
        print_keys=[]
        for i in key:
            if controller.dirs[i]!='N.A.':
                keys.append(i)
                print_keys.append('Day '+str(i))
        bsifile= np.load(controller.filepath+'/qeeg/bsi.npy')
        
        mean_alpha_a_right_pre=np.zeros(5)
        mean_alpha_a_right_post=np.zeros(5)
        mean_alpha_a_left_post=np.zeros(5)
        mean_alpha_a_left_pre=np.zeros(5)
        mean_beta_a_right_pre=np.zeros(5)
        mean_beta_a_right_post=np.zeros(5)
        mean_beta_a_left_post=np.zeros(5)
        mean_beta_a_left_pre=np.zeros(5)
        mean_delta_a_right_pre=np.zeros(5)
        mean_delta_a_right_post=np.zeros(5)
        mean_delta_a_left_post=np.zeros(5)
        mean_delta_a_left_pre=np.zeros(5)
        mean_theta_a_right_pre=np.zeros(5) 
        mean_theta_a_right_post=np.zeros(5)
        mean_theta_a_left_post=np.zeros(5)
        mean_theta_a_left_pre=np.zeros(5)
        
        (marpre,marpost,malpost,malpre,mbrpre,mbrpost,mblpost,mblpre,        mdrpre,mdrpost,mdlpost,mdlpre,mtrpre,mtrpost,mtlpost,mtlpre)=self.compute_ratios(controller, ab)
        activity='a'
        if ab==1:
            activity='b'
        j=0
        avg_bsi=np.zeros(10)
        for i in range(0, len(keys)):
            if keys[i]==1:
                j=0
            elif keys[i]==2:
                j=1
            elif keys[i]==7:
                j=2
            elif keys[i]==30:
                j=3
            else:
                j=4
            mean_alpha_a_right_pre[j]=marpre[i]
            mean_alpha_a_right_post[j]=marpost[i]
            mean_alpha_a_left_post[j]=malpost[i]
            mean_alpha_a_left_pre[j]=malpre[i]
            mean_beta_a_right_pre[j]=mbrpre[i]
            mean_beta_a_right_post[j]=mbrpost[i]
            mean_beta_a_left_post[j]=mblpost[i]
            mean_beta_a_left_pre[j]=mblpre[i]
            mean_delta_a_right_pre[j]=mdrpre[i]
            mean_delta_a_right_post[j]=mdrpost[i]
            mean_delta_a_left_post[j]=mdlpost[i]
            mean_delta_a_left_pre[j]=mdlpre[i]
            mean_theta_a_right_pre[j]=mtrpre[i]
            mean_theta_a_right_post[j]=mtrpost[i]
            mean_theta_a_left_post[j]=mtlpost[i]
            mean_theta_a_left_pre[j]=mtlpre[i]
            
            ##bsi
            plt.rcParams["figure.figsize"] = [25,16]
            plt.rcParams.update({'font.size': 12})
            bsi= bsifile[i]
            avg_bsi[2*j]= np.average(np.array( bsi[2*ab]))
            avg_bsi[2*j+1]= np.average(np.array( bsi[2*ab+1]))
            plt.plot(np.arange(len(bsi[2*ab]))*512.0/256,np.array( bsi[2*ab]),'-r', label='pre exercise')#2*1=2(b),2*0=0(a)
            plt.plot(np.arange(len(bsi[2*ab+1]))*512.0/256,np.array( bsi[2*ab+1]),'-b', label='post exercise')
            plt.xlabel('Time')
            plt.ylabel('BSI')
            plt.title('Brain Symmetry Index')
            plt.legend(loc='upper right')#2*1+1=3(bb),2*0+1=0(aa)
            plt.savefig(controller.filepath+'/qeeg/bsi_DAY_'+str(keys[i])+activity+'_.png')
            plt.close()
            
        for i in range(0, 5):
            dayframe=self.dayframes[i]
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_delta_a_right_pre[i], font=number_font)
            inFileLbl.grid(row=3, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_theta_a_right_pre[i],font=number_font)
            inFileLbl.grid(row=6, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_alpha_a_right_pre[i],font=number_font)
            inFileLbl.grid(row=10, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_beta_a_right_pre[i],font=number_font)
            inFileLbl.grid(row=14, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.2f' % avg_bsi[2*i] ,font=number_font)
            inFileLbl.grid(row=50, column=2, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.2f' % avg_bsi[2*i+1] ,font=number_font)
            inFileLbl.grid(row=50, column=4, sticky='WE', padx=0, pady=2)
            
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_delta_a_right_post[i], font=number_font)
            inFileLbl.grid(row=3, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe,text='%.1f' % mean_theta_a_right_post[i],font=number_font)
            inFileLbl.grid(row=6, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_alpha_a_right_post[i],font=number_font)
            inFileLbl.grid(row=10, column=3, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_beta_a_right_post[i],font=number_font)
            inFileLbl.grid(row=14, column=3, sticky='WE', padx=0, pady=2)
            
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_delta_a_left_pre[i], font=number_font)
            inFileLbl.grid(row=3, column=4, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_theta_a_left_pre[i],font=number_font)
            inFileLbl.grid(row=6, column=4, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_alpha_a_left_pre[i],font=number_font)
            inFileLbl.grid(row=10, column=4, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_beta_a_left_pre[i],font=number_font)
            inFileLbl.grid(row=14, column=4, sticky='WE', padx=0, pady=2)
            
            inFileLbl = tk.Label(dayframe,  text='%.1f' % mean_delta_a_left_post[i], font=number_font)
            inFileLbl.grid(row=3, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_theta_a_left_post[i],font=number_font)
            inFileLbl.grid(row=6, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_alpha_a_left_post[i],font=number_font)
            inFileLbl.grid(row=10, column=5, sticky='WE', padx=0, pady=2)
            inFileLbl = tk.Label(dayframe, text='%.1f' % mean_beta_a_left_post[i],font=number_font)
            inFileLbl.grid(row=14, column=5, sticky='WE', padx=0, pady=2)
        
        keysm=[1, 2, 3, 4, 5]
        plt.rcParams["figure.figsize"] = [25,16]
        plt.rcParams.update({'font.size': 20})
        plt.xticks(keysm, print_keys)
        #plt.bar(keysm, np.array(malpost)*np.array(marpre)/(np.array(malpre)*np.array(marpost)))
        plt.plot(keysm, np.array(malpost)*np.array(marpre)/(np.array(malpre)*np.array(marpost)))
        plt.ylabel('Alpha Ratio post vs pre left vs right')
        plt.xlabel('Day')
        plt.title('Alpha Ratio')
        plt.savefig(controller.filepath+'/qeeg/Alpha_'+activity+'_left%right_post%pre.png')
        plt.close()
        plt.xticks(keysm, print_keys)
        #plt.bar(keysm, np.array(mblpost)*np.array(mbrpre)/(np.array(mblpre)*np.array(mbrpost)))
        plt.plot(keysm, np.array(mblpost)*np.array(mbrpre)/(np.array(mblpre)*np.array(mbrpost)))
        plt.ylabel('Beta Ratio post vs pre left vs right')
        plt.xlabel('Day')
        plt.title('Beta Ratio')
        plt.savefig(controller.filepath+'/qeeg/Beta_'+activity+'_left%right_post%pre.png')
        plt.close()
        plt.xticks(keysm, print_keys)
        #plt.bar(keysm, np.array(mdlpost)*np.array(mdrpre)/(np.array(mdlpre)*np.array(mdrpost)))
        plt.plot(keysm, np.array(mdlpost)*np.array(mdrpre)/(np.array(mdlpre)*np.array(mdrpost)))
        plt.ylabel('Delta Ratio post vs pre left vs right')
        plt.xlabel('Day')
        plt.title('Delta Ratio')
        plt.savefig(controller.filepath+'/qeeg/Delta_'+activity+'_left%right_post%pre.png')
        plt.close()
        plt.xticks(keysm, print_keys)
        #plt.bar(keysm, np.array(mtlpost)*np.array(mtrpre)/(np.array(mtlpre)*np.array(mtrpost)))
        plt.plot(keysm, np.array(mtlpost)*np.array(mtrpre)/(np.array(mtlpre)*np.array(mtrpost)))
        plt.savefig(controller.filepath+'/qeeg/Theta_'+activity+'_left%right_post%pre.png')
        plt.ylabel('Theta Ratio post vs pre left vs right')
        plt.xlabel('Day')
        plt.title('Theta Ratio')
        plt.close()
        plt.xticks(keysm, print_keys)
        plt.plot(keysm, np.array(mdlpost)*np.array(marpost)/(np.array(mdrpost)*np.array(malpost)),'-b', label='post exercise')
        plt.plot(keysm, np.array(mdlpre)*np.array(marpre)/(np.array(mdrpre)*np.array(malpre)),'-r', label='pre exercise')
        plt.legend(loc='upper right')
        plt.title('Delta/Alpha Ratio')
        plt.ylabel('Delta/Alpha ratio Left vs Right')
        plt.xlabel('Day')
        plt.savefig(controller.filepath+'/qeeg/Delta%Alpha_'+activity+'_left%right_post%pre.png')
        plt.close()
        
        buttond= ttk.Button(self, text='Delta\nleft/right\npost/pre exercise',
                           command=lambda: self.show_graph(controller, activity, 'delta'))
        buttond.grid(row=70, column=15,sticky='NS', padx=2, pady=5)
        buttont= ttk.Button(self, text='Theta\nleft/right\npost/pre exercise',
                           command=lambda: self.show_graph(controller, activity, 'theta'))
        buttont.grid(row=70, column=16,sticky='NS', padx=2, pady=5)
        buttona= ttk.Button(self, text='Alpha\nleft/right\npost/pre exercise',
                           command=lambda: self.show_graph(controller, activity, 'alpha'))
        buttona.grid(row=70, column=18,sticky='NS', padx=3, pady=5)
        buttonb= ttk.Button(self, text='Beta\nleft/right\npost/pre exercise',
                           command=lambda: self.show_graph(controller, activity, 'beta'))
        buttonb.grid(row=70, column=19,sticky='NS', padx=2, pady=5)
        buttonda= ttk.Button(self, text='Delta/Alpha\nleft/right\nRatio',
                           command=lambda: self.show_graph(controller, activity, 'DeltaAlpha'))
        buttonda.grid(row=70, column=21,sticky='NS', padx=2, pady=5)
        cols= [15, 16, 18, 19, 21]

        button= ttk.Button(self, text='BSI_DAY_1',
                       command=lambda: self.show_graph(controller, activity, 'DAY_1'))
        button.grid(row=71, column=15,sticky='NS', padx=2, pady=5)
        button= ttk.Button(self, text='BSI_DAY_2',
                       command=lambda: self.show_graph(controller, activity, 'DAY_2'))
        button.grid(row=71, column=16,sticky='NS', padx=2, pady=5)
        button= ttk.Button(self, text='BSI_DAY_7',
                       command=lambda: self.show_graph(controller, activity, 'DAY_7'))
        button.grid(row=71, column=18,sticky='NS', padx=2, pady=5)
        button= ttk.Button(self, text='BSI_DAY_30',
                       command=lambda: self.show_graph(controller, activity, 'DAY_30'))
        button.grid(row=71, column=19,sticky='NS', padx=2, pady=5)
        button= ttk.Button(self, text='BSI_DAY_90',
                       command=lambda: self.show_graph(controller, activity, 'DAY_90'))
        button.grid(row=71, column=21,sticky='NS', padx=2, pady=5)
            
            
    def show_graph(self, controller, activity, band):
        if band=='delta':
            img= Image.open(controller.filepath+'/qeeg/Delta_'+activity+'_left%right_post%pre.png')
        elif band=='theta':
            img= Image.open(controller.filepath+'/qeeg/Theta_'+activity+'_left%right_post%pre.png')
        elif band=='alpha':
            img= Image.open(controller.filepath+'/qeeg/Alpha_'+activity+'_left%right_post%pre.png')
        elif band=='beta':
            img= Image.open(controller.filepath+'/qeeg/Beta_'+activity+'_left%right_post%pre.png')
        elif band=='DeltaAlpha':
            img= Image.open(controller.filepath+'/qeeg/Delta%Alpha_'+activity+'_left%right_post%pre.png')
            
        elif band== 'DAY_1':
            img= Image.open(controller.filepath+'/qeeg/bsi_DAY_1'+activity+'_.png')
            
        elif band== 'DAY_2':
            img= Image.open(controller.filepath+'/qeeg/bsi_DAY_2'+activity+'_.png')
            
        elif band== 'DAY_7':
            img= Image.open(controller.filepath+'/qeeg/bsi_DAY_7'+activity+'_.png')
        elif band== 'DAY_30':
            img= Image.open(controller.filepath+'/qeeg/bsi_DAY_30'+activity+'_.png')
        else:
            img= Image.open(controller.filepath+'/qeeg/bsi_DAY_90'+activity+'_.png')
        img.show()
        
        
    def compute_ratios(self, controller, ab):
        path=controller.filepath
        all_power= np.load(path+'/qeeg/all_powers.npy')
        spr_channels=[0,2,7,9,11,13,15,62,25,27,29,31,33,63,43,45,47,49,51,57,59]
        spr_channels=np.array(spr_channels)
        right_hem=[2, 13, 15, 31, 33, 63, 49, 51, 59]
        left_hem=[0, 7, 9, 62, 25, 27, 43, 45, 57]

        mean_right_hem=np.zeros(np.shape(all_power[0]))
        mean_left_hem=np.zeros(np.shape(all_power[0]))
        for i in right_hem:
            r=np.where(spr_channels==i)
            mean_right_hem=mean_right_hem+all_power[r[0][0]]/len(right_hem)

        for i in left_hem:
            l=np.where(spr_channels==i)
            mean_left_hem=mean_left_hem+all_power[l[0][0]]/len(left_hem)
        
        
        

        mean_alpha_a_right_pre=[] 
        mean_alpha_a_right_post=[]
        mean_alpha_a_left_post=[]
        mean_alpha_a_left_pre=[]
        mean_beta_a_right_pre=[] 
        mean_beta_a_right_post=[]
        mean_beta_a_left_post=[]
        mean_beta_a_left_pre=[]
        mean_delta_a_right_pre=[] 
        mean_delta_a_right_post=[]
        mean_delta_a_left_post=[]
        mean_delta_a_left_pre=[]
        mean_theta_a_right_pre=[] 
        mean_theta_a_right_post=[]
        mean_theta_a_left_post=[]
        mean_theta_a_left_pre=[]
        for i in range(0, len(mean_right_hem[0]), 2):

            mean_alpha_a_right_pre.append(mean_right_hem[ab][i][2])
            mean_alpha_a_left_pre.append(mean_left_hem[ab][i][2])
            mean_alpha_a_left_post.append(mean_left_hem[ab][i+1][2])
            mean_alpha_a_right_post.append(mean_right_hem[ab][i+1][2])

            mean_beta_a_right_pre.append(mean_right_hem[ab][i][3])
            mean_beta_a_left_pre.append(mean_left_hem[ab][i][3])
            mean_beta_a_left_post.append(mean_left_hem[ab][i+1][3])
            mean_beta_a_right_post.append(mean_right_hem[ab][i+1][3])

            mean_delta_a_right_pre.append(mean_right_hem[ab][i][0])
            mean_delta_a_left_pre.append(mean_left_hem[ab][i][0])
            mean_delta_a_left_post.append(mean_left_hem[ab][i+1][0])
            mean_delta_a_right_post.append(mean_right_hem[ab][i+1][0])

            mean_theta_a_right_pre.append(mean_right_hem[ab][i][1])
            mean_theta_a_left_pre.append(mean_left_hem[ab][i][1])
            mean_theta_a_left_post.append(mean_left_hem[ab][i+1][1])
            mean_theta_a_right_post.append(mean_right_hem[ab][i+1][1])

        return (mean_alpha_a_right_pre,
        mean_alpha_a_right_post,
        mean_alpha_a_left_post,
        mean_alpha_a_left_pre,
        mean_beta_a_right_pre,
        mean_beta_a_right_post,
        mean_beta_a_left_post,
        mean_beta_a_left_pre,
        mean_delta_a_right_pre ,
        mean_delta_a_right_post,
        mean_delta_a_left_post,
        mean_delta_a_left_pre,
        mean_theta_a_right_pre ,
        mean_theta_a_right_post,
        mean_theta_a_left_post,
        mean_theta_a_left_pre)


# In[20]:


app= Qeeg()
app.mainloop()


# In[6]:


print '%.1f'%3.4566


# In[7]:


#            inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_delta_a_right_pre[i])), font=number_font)
 #           inFileLbl.grid(row=3, column=2, sticky='WE', padx=0, pady=2)
  #          inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_theta_a_right_pre[i])),font=number_font)
   #         inFileLbl.grid(row=6, column=2, sticky='WE', padx=0, pady=2)
    #        inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_alpha_a_right_pre[i])),font=number_font)
     #       inFileLbl.grid(row=10, column=2, sticky='WE', padx=0, pady=2)
      #      inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_beta_a_right_pre[i])),font=number_font)
       #     inFileLbl.grid(row=14, column=2, sticky='WE', padx=0, pady=2)
            
        #    inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_delta_a_right_post[i])), font=number_font)
         #   inFileLbl.grid(row=3, column=3, sticky='WE', padx=0, pady=2)
          #  inFileLbl = tk.Label(dayframe,text='%.1f' % Decimal(str(mean_theta_a_right_post[i])),font=number_font)
           # inFileLbl.grid(row=6, column=3, sticky='WE', padx=0, pady=2)
#            inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_alpha_a_right_post[i])),font=number_font)
 #           inFileLbl.grid(row=10, column=3, sticky='WE', padx=0, pady=2)
  #          inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_beta_a_right_post[i])),font=number_font)
   #         inFileLbl.grid(row=14, column=3, sticky='WE', padx=0, pady=2)
            
    #        inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_delta_a_left_pre[i])), font=number_font)
     #       inFileLbl.grid(row=3, column=4, sticky='WE', padx=0, pady=2)
      #      inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_theta_a_left_pre[i])),font=number_font)
       #     inFileLbl.grid(row=6, column=4, sticky='WE', padx=0, pady=2)
        #    inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_alpha_a_left_pre[i])),font=number_font)
         #   inFileLbl.grid(row=10, column=4, sticky='WE', padx=0, pady=2)
          #  inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_beta_a_left_pre[i])),font=number_font)
           # inFileLbl.grid(row=14, column=4, sticky='WE', padx=0, pady=2)
            
#            inFileLbl = tk.Label(dayframe,  text='%.1f' % Decimal(str(mean_delta_a_left_post[i])), font=number_font)
 #           inFileLbl.grid(row=3, column=5, sticky='WE', padx=0, pady=2)
  #          inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_theta_a_left_post[i])),font=number_font)
   #         inFileLbl.grid(row=6, column=5, sticky='WE', padx=0, pady=2)
    #        inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_alpha_a_left_post[i])),font=number_font)
     #       inFileLbl.grid(row=10, column=5, sticky='WE', padx=0, pady=2)
      #      inFileLbl = tk.Label(dayframe, text='%.1f' % Decimal(str(mean_beta_a_left_post[i])),font=number_font)
       #     inFileLbl.grid(row=14, column=5, sticky='WE', padx=0, pady=2)

