#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 11:52:56 2024

@author: yqb22156
"""

import numpy as np
import matplotlib.pyplot as plt

in_c = 30
stck = [1.0,0.85,0.7,0.55]
crst = [1.0,1.316,1.632,1.948,2.264]
N_stck = len(stck)
N_crst = len(crst)
reps = 3

hist_crst = [[[] for _ in range(N_crst)] for _ in range(N_stck)]
hist_stck = [[[] for _ in range(N_crst)] for _ in range(N_stck)]

max_crst = 30

#average stacking

ave_stck_prop = np.zeros((N_stck,N_crst),dtype = float)
ave_crst_prop = np.zeros((N_stck,N_crst),dtype = float)
min_en_stck = np.zeros((N_stck,N_crst),dtype = float)
min_en_crst = np.zeros((N_stck,N_crst),dtype = float)

for l in range(N_stck) :
    for m in range(N_crst) :
        min_en_stck[l,m] = -1.*stck[l]*1.40726
        min_en_crst[l,m] = -1.*crst[m]*0.2375
        
tot_min_en = min_en_stck+min_en_crst

counts1 = 0

for l in range(N_stck) :
    for m in range(N_crst) :
        for k in range(reps) :
            path="STCK"+str(l)+"/CRST"+str(m)+"/"+str(k)
            
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
                
                hist_stck[l][m].append(float(vals[2])/float(vals[1]))
                ave_stck_prop[l,m] += float(vals[2])/float(vals[1])
                
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
                
                ave_crst_prop[l,m] += float(vals[2])/float(max_crst)
                
                hist_crst[l][m].append(float(vals[2])/float(max_crst))
                
            ifile.close()
            
for l in range(N_stck) :
    for m in range(N_crst) :
        ave_stck_prop[l,m] /= counts1
        ave_crst_prop[l,m] /= counts1      
   
precis = 4
  
####################### Plot Stacking Propensity

fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([0., 0., 1.0, 1.0])

plt.title(r"Stacking Propensity (24 C)")
ax.title.set_fontsize(15)

ax.set_xlabel(r"stacking strength / ox2 stacking")
ax.set_ylabel(r"cross stacking strength / ox2 cross stacking")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)

for y in range(N_crst):
    for x in range(N_stck):
        s = str(round(ave_stck_prop[x, y], precis))  # added s for plt.txt and reverse x and y for data addressing
        ax.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )
        

ticks = [0,1,2,3,4]
ax.set_yticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)


ticks = [0,1,2,3]
ax.set_xticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)

im = ax.imshow(ave_stck_prop.transpose(), cmap='jet', interpolation='nearest')

plt.savefig("Stacking_propensity_cmap.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)


####################### Plot Cross Stacking Propensity

fig = plt.figure(figsize=(5,5))
ax =fig.add_axes([0., 0., 1.0, 1.0])

plt.title(r"Cross Stacking Propensity (24 C)")
ax.title.set_fontsize(15)

ax.set_xlabel(r"stacking strength / ox2 stacking")
ax.set_ylabel(r"cross stacking strength / ox2 cross stacking")

ax.yaxis.label.set_fontsize(15)
ax.xaxis.label.set_fontsize(15)

for y in range(N_crst):
    for x in range(N_stck):
        s = str(round(ave_crst_prop[x, y], precis))  # added s for plt.txt and reverse x and y for data addressing
        plt.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )
        
ticks = [0,1,2,3,4]
ax.set_yticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)


ticks = [0,1,2,3]
ax.set_xticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)
        
im = ax.imshow(ave_crst_prop.transpose(), cmap='jet', interpolation='nearest')

plt.savefig("Cross_stacking_propensity_cmap.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)

####################### Plot Stacking Min Energy

fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
plt.title(r"Min Stacking Energy (24 C)")
ax.title.set_fontsize(15)


ax.set_ylabel(r"stacking strength / ox2 stacking")
ax.set_xlabel(r"cross stacking strength / ox2 cross stacking")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)

for y in range(N_stck):
    for x in range(N_crst):
        s = str(round(min_en_stck[y, x], precis))  # added s for plt.txt and reverse x and y for data addressing
        plt.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )
        
ticks = [0,1,2,3,4]
ax.set_xticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)


ticks = [0,1,2,3]
ax.set_yticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)
        
im = ax.imshow(min_en_stck, cmap='jet', interpolation='nearest')

plt.savefig("Stacking_minEnergy_cmap.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)

####################### Plot Cross Stacking Min Energy

fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
plt.title(r"Min Cross Stacking Energy")
ax.title.set_fontsize(15)


ax.set_ylabel(r"stacking strength / ox2 stacking")
ax.set_xlabel(r"cross stacking strength / ox2 cross stacking")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)

for y in range(N_stck):
    for x in range(N_crst):
        s = str(round(min_en_crst[y, x], precis))  # added s for plt.txt and reverse x and y for data addressing
        plt.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )
        
ticks = [0,1,2,3,4]
ax.set_xticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)


ticks = [0,1,2,3]
ax.set_yticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)
        
im = ax.imshow(min_en_crst, cmap='jet', interpolation='nearest')

plt.savefig("Cross_stacking_minEnergy_cmap.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)

####################### Plot Cross Stacking Min Energy

fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
plt.title(r"Min Stacking + Min Cross Stacking Energy (24 C)")
ax.title.set_fontsize(15)


ax.set_xlabel(r"stacking strength / ox2 stacking")
ax.set_ylabel(r"cross stacking strength / ox2 cross stacking")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)

for y in range(N_crst):
    for x in range(N_stck):
        s = str(round(tot_min_en[x, y], precis))  # added s for plt.txt and reverse x and y for data addressing
        plt.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )
        
ticks = [0,1,2,3,4]
ax.set_yticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)


ticks = [0,1,2,3]
ax.set_xticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)
        
im = ax.imshow(tot_min_en.transpose(), cmap='jet', interpolation='nearest')

plt.savefig("Sum_minEnergy_cmap.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)



#histograms


fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
plt.title(r"Histogram of Cr Stacked pairs (24 C)")
ax.title.set_fontsize(15)


ax.set_ylabel(r"probability density")
ax.set_xlabel(r"fraction of cross stacked pairs")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)
"""
for y in range(N_stck):
    for x in range(N_crst):
        s = str(round(tot_min_en[y, x], precis))  # added s for plt.txt and reverse x and y for data addressing
        plt.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )

        
ticks = [0,1,2,3,4]
ax.set_xticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)


ticks = [0,1,2,3]
ax.set_yticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)
"""       

"""
for l in range(N_stck) :
    for m in range(N_crst) :
        for k in range(len(hist_crst[l][m])) :
            hist_crst[l][m][k] /= len(hist_crst[l][m])
"""

bins = []

for k in range(max_crst):
    bins.append((-0.5+k)/max_crst)

print(bins)

#print(hist_crst[0][4])


plt.hist(hist_crst[0][4], bins, label = "stck = 1.0, crst = 2.26", density = True)
plt.hist(hist_crst[0][3], bins, label = "stck = 1.0, crst = 1.95",  density = True)
plt.hist(hist_crst[0][2], bins, label = "stck = 1.0, crst = 1.63",  density = True)
plt.hist(hist_crst[0][1], bins, label = "stck = 1.0, crst = 1.32",  density = True)
plt.hist(hist_crst[0][0], bins, label = "stck = 1.0, crst = 1.0",  density = True)



ax.legend(loc='upper right', prop={'size':15})


plt.savefig("Histo_cross_ints.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)


#2

fig = plt.figure(figsize=(5,5))

ax =fig.add_axes([0., 0., 1.0, 1.0])
plt.title(r"Histogram of Cr Stacked pairs (24 C)")
ax.title.set_fontsize(15)


ax.set_ylabel(r"probability density")
ax.set_xlabel(r"fraction of cross stacked pairs")

ax.xaxis.label.set_fontsize(15)
ax.yaxis.label.set_fontsize(15)
"""
for y in range(N_stck):
    for x in range(N_crst):
        s = str(round(tot_min_en[y, x], precis))  # added s for plt.txt and reverse x and y for data addressing
        plt.text(x, y, s,
            horizontalalignment='center',
            verticalalignment='center',
            color='white',
            )

        
ticks = [0,1,2,3,4]
ax.set_xticks(ticks)

dic = {0 : "$"+str(crst[0])+"$", 1 : "$"+str(crst[1])+"$", 2 : "$"+str(crst[2])+"$", 3 : "$"+str(crst[3])+"$",4 : "$"+str(crst[4])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_xticklabels(labels)


ticks = [0,1,2,3]
ax.set_yticks(ticks)

dic = {0 : "$"+str(stck[0])+"$", 1 : "$"+str(stck[1])+"$", 2 : "$"+str(stck[2])+"$", 3 : "$"+str(stck[3])+"$"}
labels = [ticks[i] if t not in dic.keys() else dic[t] for i,t in enumerate(ticks)]
## or 
# labels = [dic.get(t, ticks[i]) for i,t in enumerate(ticks)]
ax.set_yticklabels(labels)
"""       

"""
for l in range(N_stck) :
    for m in range(N_crst) :
        for k in range(len(hist_crst[l][m])) :
            hist_crst[l][m][k] /= len(hist_crst[l][m])
"""

bins = []

for k in range(max_crst):
    bins.append((-0.5+k)/max_crst)

print(bins)

#print(hist_crst[0][4])


plt.hist(hist_crst[1][4], bins, label = "stck = 0.85, crst = 2.26", density = True)
plt.hist(hist_crst[1][3], bins, label = "stck = 0.85, crst = 1.95",  density = True)
plt.hist(hist_crst[1][2], bins, label = "stck = 0.85, crst = 1.63",  density = True)
plt.hist(hist_crst[1][1], bins, label = "stck = 0.85, crst = 1.32",  density = True)
plt.hist(hist_crst[1][0], bins, label = "stck = 0.85, crst = 1.0",  density = True)



ax.legend(loc='upper right', prop={'size':15})


plt.savefig("Histo_cross_ints_stck085.pdf", bbox_inches = 'tight',  pad_inches = 0, transparent=True)