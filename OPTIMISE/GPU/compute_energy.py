#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:04:54 2024

@author: yqb22156
"""

import numpy as np
import copy
import math

#arrays containing relative coordinates of all interacting pairs.
#BONDED
# 0: r fene
# 1: r stck
# 2: costh4
# 3: costh5
# 4: costh6
# 5: cosphi1
# 6: cosphi2
# 7,8,9,10: types

#NOBOND
# 0: r hydr
# 1: costh1
# 2: costh2
# 3: costh4
# 4: costh7
# 5: costh8
# 6,7,8,9: types


trj_bond_pairs_coord = []
trj_nobond_pairs_coord = []

# Parameters is a 5D tensor. i,j,k,l,m 
# i: parameter id
# j,k,l,m : bases



PARS_LIST = [
            'FENE_EPS', 'FENE_R0', 'FENE_DELTA', 'FENE_DELTA2', 
             'HYDR_EPS', 'HYDR_R0', 'HYDR_A', 'HYDR_RC', 'HYDR_BLOW', 'HYDR_BHIGH', 'HYDR_RLOW', 'HYDR_RHIGH', 'HYDR_RCLOW', 'HYDR_RCHIGH',
             'HYDR_THETA1_T0', 'HYDR_THETA1_A', 'HYDR_THETA1_B', 'HYDR_THETA1_TS', 'HYDR_THETA1_TC',
             'HYDR_THETA2_T0', 'HYDR_THETA2_A', 'HYDR_THETA2_B', 'HYDR_THETA2_TS', 'HYDR_THETA2_TC',
             'HYDR_THETA3_T0', 'HYDR_THETA3_A', 'HYDR_THETA3_B', 'HYDR_THETA3_TS', 'HYDR_THETA3_TC',
             'HYDR_THETA4_T0', 'HYDR_THETA4_A', 'HYDR_THETA4_B', 'HYDR_THETA4_TS', 'HYDR_THETA4_TC',
             'HYDR_THETA7_T0', 'HYDR_THETA7_A', 'HYDR_THETA7_B', 'HYDR_THETA7_TS', 'HYDR_THETA7_TC',
             'HYDR_THETA8_T0', 'HYDR_THETA8_A', 'HYDR_THETA8_B', 'HYDR_THETA8_TS', 'HYDR_THETA8_TC',
             'STCK_EPS', 'STCK_R0', 'STCK_A', 'STCK_RC', 'STCK_BLOW', 'STCK_BHIGH', 'STCK_RLOW', 'STCK_RHIGH', 'STCK_RCLOW', 'STCK_RCHIGH',
             'STCK_THETA4_T0', 'STCK_THETA4_A', 'STCK_THETA4_B', 'STCK_THETA4_TS', 'STCK_THETA4_TC',
             'STCK_THETA5_T0', 'STCK_THETA5_A', 'STCK_THETA5_B', 'STCK_THETA5_TS', 'STCK_THETA5_TC',
             'STCK_THETA6_T0', 'STCK_THETA6_A', 'STCK_THETA6_B', 'STCK_THETA6_TS', 'STCK_THETA6_TC',
             'STCK_PHI1_A', 'STCK_PHI1_B', 'STCK_PHI1_XC', 'STCK_PHI1_XS',
             'STCK_PHI2_A', 'STCK_PHI2_B', 'STCK_PHI2_XC', 'STCK_PHI2_XS',
             'CRST_K_33', 'CRST_R0_33', 'CRST_RC_33', 'CRST_BLOW_33', 'CRST_BHIGH_33', 'CRST_RLOW_33', 'CRST_RHIGH_33', 'CRST_RCLOW_33', 'CRST_RCHIGH_33',
             'CRST_THETA1_T0_33', 'CRST_THETA1_A_33', 'CRST_THETA1_B_33', 'CRST_THETA1_TS_33', 'CRST_THETA1_TC_33',
             'CRST_THETA2_T0_33', 'CRST_THETA2_A_33', 'CRST_THETA2_B_33', 'CRST_THETA2_TS_33', 'CRST_THETA2_TC_33',
             'CRST_THETA3_T0_33', 'CRST_THETA3_A_33', 'CRST_THETA3_B_33', 'CRST_THETA3_TS_33', 'CRST_THETA3_TC_33',
             'CRST_THETA4_T0_33', 'CRST_THETA4_A_33', 'CRST_THETA4_B_33', 'CRST_THETA4_TS_33', 'CRST_THETA4_TC_33',
             'CRST_THETA7_T0_33', 'CRST_THETA7_A_33', 'CRST_THETA7_B_33', 'CRST_THETA7_TS_33', 'CRST_THETA7_TC_33',
             'CRST_THETA8_T0_33', 'CRST_THETA8_A_33', 'CRST_THETA8_B_33', 'CRST_THETA8_TS_33', 'CRST_THETA8_TC_33',
             'CRST_K_55', 'CRST_R0_55', 'CRST_RC_55', 'CRST_BLOW_55', 'CRST_BHIGH_55', 'CRST_RLOW_55', 'CRST_RHIGH_55', 'CRST_RCLOW_55', 'CRST_RCHIGH_55',
             'CRST_THETA1_T0_55', 'CRST_THETA1_A_55', 'CRST_THETA1_B_55', 'CRST_THETA1_TS_55', 'CRST_THETA1_TC_55',
             'CRST_THETA2_T0_55', 'CRST_THETA2_A_55', 'CRST_THETA2_B_55', 'CRST_THETA2_TS_55', 'CRST_THETA2_TC_55',
             'CRST_THETA3_T0_55', 'CRST_THETA3_A_55', 'CRST_THETA3_B_55', 'CRST_THETA3_TS_55', 'CRST_THETA3_TC_55',
             'CRST_THETA4_T0_55', 'CRST_THETA4_A_55', 'CRST_THETA4_B_55', 'CRST_THETA4_TS_55', 'CRST_THETA4_TC_55',
             'CRST_THETA7_T0_55', 'CRST_THETA7_A_55', 'CRST_THETA7_B_55', 'CRST_THETA7_TS_55', 'CRST_THETA7_TC_55',
             'CRST_THETA8_T0_55', 'CRST_THETA8_A_55', 'CRST_THETA8_B_55', 'CRST_THETA8_TS_55', 'CRST_THETA8_TC_55'
             ]



#create list with indices of parameters used in compute_energy.
#The purpuse of this is that even if PARS_LIST gets changed, compute_energy doesn't need to be updated. This should be faster than searching within a map
#ids are reported below
par_index = []

    
#FENE
    
par_index.append(PARS_LIST.index("FENE_EPS"))
par_index.append(PARS_LIST.index("FENE_R0"))
par_index.append(PARS_LIST.index("FENE_DELTA"))
par_index.append(PARS_LIST.index("FENE_DELTA2"))


#HYDROGEN

par_index.append(PARS_LIST.index("HYDR_EPS"))

par_index.append(PARS_LIST.index("HYDR_R0"))
par_index.append(PARS_LIST.index("HYDR_A"))
par_index.append(PARS_LIST.index("HYDR_RC"))
par_index.append(PARS_LIST.index("HYDR_BLOW"))
par_index.append(PARS_LIST.index("HYDR_BHIGH"))
par_index.append(PARS_LIST.index("HYDR_RLOW"))
par_index.append(PARS_LIST.index("HYDR_RHIGH"))
par_index.append(PARS_LIST.index("HYDR_RCLOW"))
par_index.append(PARS_LIST.index("HYDR_RCHIGH"))

par_index.append(PARS_LIST.index("HYDR_THETA1_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA1_A"))
par_index.append(PARS_LIST.index("HYDR_THETA1_B"))
par_index.append(PARS_LIST.index("HYDR_THETA1_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA1_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA2_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA2_A"))
par_index.append(PARS_LIST.index("HYDR_THETA2_B"))
par_index.append(PARS_LIST.index("HYDR_THETA2_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA2_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA3_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA3_A"))
par_index.append(PARS_LIST.index("HYDR_THETA3_B"))
par_index.append(PARS_LIST.index("HYDR_THETA3_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA3_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA4_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA4_A"))
par_index.append(PARS_LIST.index("HYDR_THETA4_B"))
par_index.append(PARS_LIST.index("HYDR_THETA4_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA4_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA4_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA4_A"))
par_index.append(PARS_LIST.index("HYDR_THETA4_B"))
par_index.append(PARS_LIST.index("HYDR_THETA4_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA4_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA7_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA7_A"))
par_index.append(PARS_LIST.index("HYDR_THETA7_B"))
par_index.append(PARS_LIST.index("HYDR_THETA7_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA7_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA8_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA8_A"))
par_index.append(PARS_LIST.index("HYDR_THETA8_B"))
par_index.append(PARS_LIST.index("HYDR_THETA8_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA8_TC"))


#STACKING

par_index.append(PARS_LIST.index("STCK_EPS"))

par_index.append(PARS_LIST.index("STCK_R0"))
par_index.append(PARS_LIST.index("STCK_A"))
par_index.append(PARS_LIST.index("STCK_RC"))
par_index.append(PARS_LIST.index("STCK_BLOW"))
par_index.append(PARS_LIST.index("STCK_BHIGH"))
par_index.append(PARS_LIST.index("STCK_RLOW"))
par_index.append(PARS_LIST.index("STCK_RHIGH"))
par_index.append(PARS_LIST.index("STCK_RCLOW"))
par_index.append(PARS_LIST.index("STCK_RCHIGH"))

par_index.append(PARS_LIST.index("STCK_THETA4_T0"))
par_index.append(PARS_LIST.index("STCK_THETA4_A"))
par_index.append(PARS_LIST.index("STCK_THETA4_B"))
par_index.append(PARS_LIST.index("STCK_THETA4_TS"))
par_index.append(PARS_LIST.index("STCK_THETA4_TC"))

par_index.append(PARS_LIST.index("STCK_THETA5_T0"))
par_index.append(PARS_LIST.index("STCK_THETA5_A"))
par_index.append(PARS_LIST.index("STCK_THETA5_B"))
par_index.append(PARS_LIST.index("STCK_THETA5_TS"))
par_index.append(PARS_LIST.index("STCK_THETA5_TC"))

par_index.append(PARS_LIST.index("STCK_THETA6_T0"))
par_index.append(PARS_LIST.index("STCK_THETA6_A"))
par_index.append(PARS_LIST.index("STCK_THETA6_B"))
par_index.append(PARS_LIST.index("STCK_THETA6_TS"))
par_index.append(PARS_LIST.index("STCK_THETA6_TC"))

par_index.append(PARS_LIST.index("STCK_PHI1_A"))
par_index.append(PARS_LIST.index("STCK_PHI1_B"))
par_index.append(PARS_LIST.index("STCK_PHI1_XC"))
par_index.append(PARS_LIST.index("STCK_PHI1_XS"))

par_index.append(PARS_LIST.index("STCK_PHI2_A"))
par_index.append(PARS_LIST.index("STCK_PHI2_B"))
par_index.append(PARS_LIST.index("STCK_PHI2_XC"))
par_index.append(PARS_LIST.index("STCK_PHI2_XS"))


#CROSS STACKING
    
par_index.append(PARS_LIST.index("CRST_K_33"))
par_index.append(PARS_LIST.index("CRST_R0_33"))
par_index.append(PARS_LIST.index("CRST_RC_33"))
par_index.append(PARS_LIST.index("CRST_BLOW_33"))
par_index.append(PARS_LIST.index("CRST_BHIGH_33"))
par_index.append(PARS_LIST.index("CRST_RLOW_33"))
par_index.append(PARS_LIST.index("CRST_RHIGH_33"))
par_index.append(PARS_LIST.index("CRST_RCLOW_33"))
par_index.append(PARS_LIST.index("CRST_RCHIGH_33"))

par_index.append(PARS_LIST.index("CRST_THETA1_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA2_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA3_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA4_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA7_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA8_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_TC_33"))

par_index.append(PARS_LIST.index("CRST_K_55"))
par_index.append(PARS_LIST.index("CRST_R0_55"))
par_index.append(PARS_LIST.index("CRST_RC_55"))
par_index.append(PARS_LIST.index("CRST_BLOW_55"))
par_index.append(PARS_LIST.index("CRST_BHIGH_55"))
par_index.append(PARS_LIST.index("CRST_RLOW_55"))
par_index.append(PARS_LIST.index("CRST_RHIGH_55"))
par_index.append(PARS_LIST.index("CRST_RCLOW_55"))
par_index.append(PARS_LIST.index("CRST_RCHIGH_55"))

par_index.append(PARS_LIST.index("CRST_THETA1_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA2_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA3_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA4_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA7_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA8_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_TC_55"))


"""
IDS:
    
FENE
    
0 - FENE_EPS
1 - FENE_R0
2 - FENE_DELTA
3 - FENE_DELTA2


HYDROGEN

4 - HYDR_EPS

5 - HYDR_R0
6 - HYDR_A
7 - HYDR_RC
8 - HYDR_BLOW
9 - HYDR_BHIGH
10 - HYDR_RLOW
11 - HYDR_RHIGH
12 - HYDR_RCLOW
13 - HYDR_RCHIGH

14 - HYDR_THETA1_T0
15 - HYDR_THETA1_A
16 - HYDR_THETA1_B
17 - HYDR_THETA1_TS
18 - HYDR_THETA1_TC

19 - HYDR_THETA2_T0
20 - HYDR_THETA2_A
21 - HYDR_THETA2_B
22 - HYDR_THETA2_TS
23 - HYDR_THETA2_TC

24 - HYDR_THETA3_T0
25 - HYDR_THETA3_A
26 - HYDR_THETA3_B
27 - HYDR_THETA3_TS
28 - HYDR_THETA3_TC

29 - HYDR_THETA4_T0
30 - HYDR_THETA4_A
31 - HYDR_THETA4_B
32 - HYDR_THETA4_TS
33 - HYDR_THETA4_TC

29 - HYDR_THETA4_T0
30 - HYDR_THETA4_A
31 - HYDR_THETA4_B
32 - HYDR_THETA4_TS
33 - HYDR_THETA4_TC

34 - HYDR_THETA7_T0
35 - HYDR_THETA7_A
36 - HYDR_THETA7_B
37 - HYDR_THETA7_TS
38 - HYDR_THETA7_TC

39 - HYDR_THETA8_T0
40 - HYDR_THETA8_A
41 - HYDR_THETA8_B
42 - HYDR_THETA8_TS
43 - HYDR_THETA8_TC


STACKING

44 - STCK_EPS

45 - STCK_R0
46 - STCK_A
47 - STCK_RC
48 - STCK_BLOW
49 - STCK_BHIGH
50 - STCK_RLOW
51 - STCK_RHIGH
52 - STCK_RCLOW
53 - STCK_RCHIGH

54 - STCK_THETA4_T0
55 - STCK_THETA4_A
56 - STCK_THETA4_B
57 - STCK_THETA4_TS
58 - STCK_THETA4_TC

59 - STCK_THETA5_T0
60 - STCK_THETA5_A
61 - STCK_THETA5_B
62 - STCK_THETA5_TS
63 - STCK_THETA5_TC

64 - STCK_THETA6_T0
65 - STCK_THETA6_A
66 - STCK_THETA6_B
67 - STCK_THETA6_TS
68 - STCK_THETA6_TC

69 - STCK_PHI1_A
70 - STCK_PHI1_B
71 - STCK_PHI1_XC
72 - STCK_PHI1_XS

73 - STCK_PHI2_A
74 - STCK_PHI2_B
75 - STCK_PHI2_XC
76 - STCK_PHI2_XS


CROSS STACKING
    
77 - CRST_K_33
78 - CRST_R0_33
79 - CRST_RC_33
80 - CRST_BLOW_33
81 - CRST_BHIGH_33
82 - CRST_RLOW_33
83 - CRST_RHIGH_33
84 - CRST_RCLOW_33
85 - CRST_RCHIGH_33

86 - CRST_THETA1_T0_33
87 - CRST_THETA1_A_33
88 - CRST_THETA1_B_33
89 - CRST_THETA1_TS_33
90 - CRST_THETA1_TC_33

91 - CRST_THETA2_T0_33
92 - CRST_THETA2_A_33
93 - CRST_THETA2_B_33
94 - CRST_THETA2_TS_33
95 - CRST_THETA2_TC_33

96 - CRST_THETA3_T0_33
97 - CRST_THETA3_A_33
98 - CRST_THETA3_B_33
99 - CRST_THETA3_TS_33
100 - CRST_THETA3_TC_33

101 - CRST_THETA4_T0_33
102 - CRST_THETA4_A_33
103 - CRST_THETA4_B_33
104 - CRST_THETA4_TS_33
105 - CRST_THETA4_TC_33

106 - CRST_THETA7_T0_33
107 - CRST_THETA7_A_33
108 - CRST_THETA7_B_33
109 - CRST_THETA7_TS_33
110 - CRST_THETA7_TC_33

111 - CRST_THETA8_T0_33
112 - CRST_THETA8_A_33
113 - CRST_THETA8_B_33
114 - CRST_THETA8_TS_33
115 - CRST_THETA8_TC_33

116 - CRST_K_55
117 - CRST_R0_55
118 - CRST_RC_55
119 - CRST_BLOW_55
120 - CRST_BHIGH_55
121 - CRST_RLOW_55
122 - CRST_RHIGH_55
123 - CRST_RCLOW_55
124 - CRST_RCHIGH_55

125 - CRST_THETA1_T0_55
126 - CRST_THETA1_A_55
127 - CRST_THETA1_B_55
128 - CRST_THETA1_TS_55
129 - CRST_THETA1_TC_55

130 - CRST_THETA2_T0_55
131 - CRST_THETA2_A_55
132 - CRST_THETA2_B_55
133 - CRST_THETA2_TS_55
134 - CRST_THETA2_TC_55

135 - CRST_THETA3_T0_55
136 - CRST_THETA3_A_55
137 - CRST_THETA3_B_55
138 - CRST_THETA3_TS_55
139 - CRST_THETA3_TC_55

140 - CRST_THETA4_T0_55
141 - CRST_THETA4_A_55
142 - CRST_THETA4_B_55
143 - CRST_THETA4_TS_55
144 - CRST_THETA4_TC_55

145 - CRST_THETA7_T0_55
146 - CRST_THETA7_A_55
147 - CRST_THETA7_B_55
148 - CRST_THETA7_TS_55
149 - CRST_THETA7_TC_55

150 - CRST_THETA8_T0_55
151 - CRST_THETA8_A_55
152 - CRST_THETA8_B_55
153 - CRST_THETA8_TS_55
154 - CRST_THETA8_TC_55

"""

bases = ['A','C','G','T']
#map bases to id
def base_to_id(b) :
    if b == 'A':
        return 0
    elif  b == 'C':
        return 1
    elif  b == 'G':
        return 2
    elif  b == 'T':
        return 3
    else:
        return -1


#MODULATIONS

def Vmod(theta, a, theta0) :
    f = 1-a*(theta-theta0)**2
    return f

def Vsmooth(x,b,xc):
    f = b*(xc-x)**2
    return f

def Morse(r,epsilon,r0,a) :
    f = epsilon*(1-math.exp(-(r-r0)*a))**2
    return f

def Vharmonic(x,x0):
    f = 0.5*(x0-x)**2
    return f

#pass shift as argument: no need to compute it more than once
def f1(r,shift,epsilon,r0,a,rc,bl,bh,rl,rh,rcl,rch) :
    f = 0
    if r >= rl and r <= rh :
        f = Morse(r,epsilon,r0,a) - shift
    elif r > rcl and r < rl :
        f = epsilon*Vsmooth(r,bl,rcl)
    elif r > rh and r < rch :
        f = epsilon*Vsmooth(r,bh,rch)
    else :
        f = 0
        
    return f

def f2(r,r0,rc,bl,bh,rl,rh,rcl,rch) :
    f = 0
    if r >= rl and r <= rh :
        f = Vharmonic(r,r0)-Vharmonic(rc,r0)
    elif r > rcl and r < rl :
        f = Vsmooth(r,bl,rcl)
    elif r > rh and r < rch :
        f = Vsmooth(r,bh,rch)
    else :
        f = 0
        
    return f

def f4(theta, a, b, theta0, dtheta_s, dtheta_c) :
    f = 0
    if theta >= (theta0-dtheta_s) and theta <= (theta0+dtheta_s) :
        f = Vmod(theta,a,theta0)
    elif theta > (theta0-dtheta_c) and theta < (theta0-dtheta_s) :
        f = Vsmooth(theta,b,theta0-dtheta_c)
    elif theta > (theta0+dtheta_s) and theta < (theta0+dtheta_c) :
        f = Vsmooth(theta,b,theta0+dtheta_c)
    else :
        f = 0
    
    return f


#5th base type = ghost. Used for ends, where there are no neighbours
OXPS_zero = np.zeros((156,5,5,5,5),dtype=float)
shifts = np.zeros((2,5,5,5,5),dtype=float) #0 = hydr, 1 = stck
OPTI_P = []
OPTI_P_LIST = []
OPTI_P_INDICES = np.zeros((len(OPTI_P_LIST),5),dtype=int)

for i in len(OPTI_P_LIST) :
    indices = []
    words = OPTI_P_LIST.split("_")
    string = words[0]
    for j in range(1,len(words)-4):
        string+="_"+words[j]
    index = PARS_LIST.index()
    if index < 0 or index > 156:
        print("Warning: parameter "+OPTI_P_LIST[i]+" not registered.")
        continue
    indices.append(index)
    indices.append(base_to_id(words[len(words)-4]))
    indices.append(base_to_id(words[len(words)-3]))
    indices.append(base_to_id(words[len(words)-2]))
    indices.append(base_to_id(words[len(words)-1]))
    
    OPTI_P_INDICES.append(indices)
 


#find initial values of the parameters from model.h file (if it gets updated, changes are read without modifying the code)
pars_from_modelh = []
vals_from_modelh = []

model_file = open("",'r')
####TODO
    
    
#over_indices and over_vals are indices and values of parameters to overwrite
#note: overwrite pars must include STCK_x_y
#T = temperature in oxdna units. This is needed to correctly set STCK_EPS
def init_oxpars(over_indices, over_vals,T) :
    
    stck_fact_eps = 1.
    
    for i in range(len(PARS_LIST)) :
        
        if PARS_LIST[i] == "STCK_EPS":
            index = pars_from_modelh.index("STCK_FACT_EPS")
            stck_fact_eps = vals_from_modelh[index]
            for j in range(5) :
                for k in range(5) :
                    for l in range(5):
                        for m in range(5) :
                            OXPS_zero[i][j][k][l][m] = 1.
            continue
        
        index = pars_from_modelh.index(PARS_LIST[i])
        val = vals_from_modelh[index]
        
        for j in range(5) :
            for k in range(5) :
                for l in range(5):
                    for m in range(5) :
                        OXPS_zero[i][j][k][l][m] = val
    
    #here we use the initial custom parameters, including stck_x_y
    for i in range(len(over_indices)) :
        OXPS_zero[over_indices[i][0]][over_indices[i][1]][over_indices[i][2]][over_indices[i][3]][over_indices[i][4]] = over_vals[i]
        
    #set eps and shifts
    for j in range(5) :
        for k in range(5) :
            for l in range(5):
                for m in range(5) :
                    #hydr
                    shifts[0][j][k][l][m] = Morse(OXPS_zero[par_index[7]][j][k][l][m],OXPS_zero[par_index[4]][j][k][l][m],OXPS_zero[par_index[5]][j][k][l][m],OXPS_zero[par_index[6]][j][k][l][m])
                    #stacking
                    OXPS_zero[par_index[44]][j][k][l][m] =  OXPS_zero[par_index[44]][j][k][l][m]* (1.0 - stck_fact_eps + (T * 9.0 * stck_fact_eps))
                    shifts[1][j][k][l][m] = Morse(OXPS_zero[par_index[47]][j][k][l][m],OXPS_zero[par_index[44]][j][k][l][m],OXPS_zero[par_index[45]][j][k][l][m],OXPS_zero[par_index[46]][j][k][l][m])
                    
    
    
    return
  


PARS = copy.deepcopy(OXPS_zero) #tensor with current value of the parameterset. Use this to compute energy 
#overwrite opti parameters
for i in range(len(OPTI_P_INDICES)):
    PARS[OPTI_P_INDICES[i][0]][OPTI_P_INDICES[i][1]][OPTI_P_INDICES[i][2]][OPTI_P_INDICES[i][3]][OPTI_P_INDICES[i][4]] = OPTI_P[i]
 




#arrays containing relative coordinates of all interacting pairs.
#BONDED
# 0: r fene
# 1: r stck
# 2: costh4
# 3: costh5
# 4: costh6
# 5: cosphi1
# 6: cosphi2
# 7,8,9,10: types

#NOBOND
# 0: r hydr
# 1: costh1
# 2: costh2
# 3: costh4
# 4: costh7
# 5: costh8
# 6,7,8,9: types


#given a configuration, computes the oxdna potential energy.
#only FENE, hydrogen, stacking and cross-stacking
#loop over bonded and unbonded pairs (identified when trajectory is load) and compute energy with parameterset = PARS


def compute_energy(bd_pairs, unbd_pairs, PARS) :
    
    energy = 0.
    
    
    #bonded
    
    for pair in bd_pairs :
        t1 = pair[7]
        t2 = pair[8]
        t3 = pair[9]
        t4 = pair[10]
        #fene
        #evenually, here we can handle case with _use_mbf = True
        
        if(math.abs(pair[0]-PARS[par_index[1]][t1][t2][t3][t4])>PARS[par_index[2]][t1][t2][t3][t4])-1.E-16 :
            energy += 1.E9
            return energy
        
        energy -= PARS[par_index[0]][t1][t2][t3][t4]/2.*math.log( 1.-math.sqrt( pair[0]-PARS[par_index[1]][t1][t2][t3][t4] )/PARS[par_index[2]][t1][t2][t3][t4] )
           
        #stacking
        energy += f1(pair[1],PARS[par_index[4]][t1][t2][t3][t4],PARS[par_index[5]][t1][t2][t3][t4],PARS[par_index[6]][t1][t2][t3][t4],PARS[par_index[7]][t1][t2][t3][t4],PARS[par_index[8]][t1][t2][t3][t4], \
                     PARS[par_index[9]][t1][t2][t3][t4],PARS[par_index[10]][t1][t2][t3][t4],PARS[par_index[11]][t1][t2][t3][t4],PARS[par_index[12]][t1][t2][t3][t4],PARS[par_index[13]][t1][t2][t3][t4])*\
                  f4(pair[2],PARS[par_index[14]][t1][t2][t3][t4],PARS[par_index[15]][t1][t2][t3][t4],PARS[par_index[16]][t1][t2][t3][t4],PARS[par_index[17]][t1][t2][t3][t4],PARS[par_index[18]][t1][t2][t3][t4] )*\
                  f4(-pair[3],PARS[par_index[19]][t1][t2][t3][t4],PARS[par_index[20]][t1][t2][t3][t4],PARS[par_index[21]][t1][t2][t3][t4],PARS[par_index[22]][t1][t2][t3][t4],PARS[par_index[23]][t1][t2][t3][t4] )*\
                  f4(pair[4],PARS[par_index[24]][t1][t2][t3][t4],PARS[par_index[25]][t1][t2][t3][t4],PARS[par_index[26]][t1][t2][t3][t4],PARS[par_index[27]][t1][t2][t3][t4],PARS[par_index[28]][t1][t2][t3][t4] )
            
        
    
    return
                    
    

