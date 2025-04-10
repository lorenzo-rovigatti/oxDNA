#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 11:52:56 2024

@author: yqb22156
"""

import numpy as np
import matplotlib.pyplot as plt
import "../../colours.py"

#style

plt.style.use("../../style.sty")

seqs=['AA', 'AG', 'AC', 'AT', 'CG', 'CC', 'CT', 'GG', 'GT', 'TT', 'RAND']

in_c = 30

Nseqs = len(poly_steps)

hist_crst = [[] for _ in range(Nseqs)]
hist_stck = [[] for _ in range(Nseqs)]

max_crst = 30

#average stacking

ave_stck_prop = np.zeros(Nseqs,dtype = float)
ave_crst_prop = np.zeros(Nseqs,dtype = float)

counts1 = 0

for l in range(Nseqs) :
            path=str(l)
            
            ifile = open(path+"/stacking_propensity.dat", "r")
            
            counts_tmp = 0
            
            for line in ifile.readlines() :
                if line[0] == '#' or line[0] == 'A' or line[0] == 'C' or line[0] == 'G' or line[0] == 'T' or line[0] == '\n':
                    continue
                
                counts_tmp += 1
                
                if l == 0 and m == 0 :
                    counts1+=1
                
                if counts_tmp < in_c :
                    continue
                
                vals = line.split()
                
                #print(l,m)
                
                hist_stck[l].append(float(vals[2])/float(vals[1]))
                ave_stck_prop[l] += float(vals[2])/float(vals[1])
                
            ifile.close()
            
            ifile = open(path+"/cross_stacking_propensity.dat", "r")
            
            counts_tmp = 0
            
            for line in ifile.readlines() :
                if line[0] == '#' or line[0] == 'A' or line[0] == 'C' or line[0] == 'G' or line[0] == 'T' or line[0] == '\n':
                    continue
                
                counts_tmp+=1
                
                if counts_tmp < in_c :
                    continue
                
                vals = line.split()
                
                ave_crst_prop[l] += float(vals[2])/float(max_crst)
                
                hist_crst[l].append(float(vals[2])/float(max_crst))
                
            ifile.close()
            
for l in range(Nseqs) :
    ave_stck_prop[l] /= counts1
    ave_crst_prop[l] /= counts1      
   
precis = 4
  
####################### Plot Stacking Propensity

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$sequence$",fontsize=20)
ax.set_ylabel(r"stacking propensity",fontsize=20)
#ax.set_ylim(30,60)

ticks = np.arrange(Nseqs)
ax.set_xticks(ticks)

dic = dict(zip(ticks,poly_steps))

labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)

#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(ticks,ave_stck_prop,color=CCS_G[1],label=None)
#ax.legend(fontsize = 20)

plt.savefig("stacking_propensity.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)


####################### Plot Cross Stacking Propensity

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$sequence$",fontsize=20)
ax.set_ylabel(r"cross stacking propensity",fontsize=20)
#ax.set_ylim(30,60)

ticks = np.arrange(Nseqs)
ax.set_xticks(ticks)

dic = dict(zip(ticks,poly_steps))

labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)

#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(ticks,ave_crst_prop,color=CCS_G[1],label=None)
#ax.legend(fontsize = 20)

plt.savefig("cross_stacking_propensity.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)


#histograms

#stck

fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
ax.title.set_fontsize(15)


ax.set_ylabel(r"probability density")
ax.set_xlabel(r"fraction of stacked pairs")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)


bins = []

for k in range(max_crst):
    bins.append((-0.5+k)/max_crst)

print(bins)

#print(hist_crst[0][4])


for i in range(Nseqs) :
    plt.hist(hist_stck[i], bins, color=CCS_G[i], label = seqs[i], density = True)


ax.legend(loc='upper right', prop={'size':15})

#crst

plt.savefig("histo_stck.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)

fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
#plt.title(r"Histogram of Cr Stacked pairs (24 C)")
ax.title.set_fontsize(15)


ax.set_ylabel(r"probability density")
ax.set_xlabel(r"fraction of cross stacked pairs")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)

bins = []

for k in range(max_crst):
    bins.append((-0.5+k)/max_crst)

print(bins)

#print(hist_crst[0][4])

for i in range(Nseqs) :
    plt.hist(hist_crst[i], bins, color=CCS_G[i], label = seqs[i], density = True)

ax.legend(loc='upper right', prop={'size':15})


plt.savefig("histo_crst.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)

