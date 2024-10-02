
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:04:54 2024

@author: yqb22156
"""

import numpy as np
import copy
import math
import torch

#list of oxdna parameters used in energy computation
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


OXPS_zero = np.zeros((len(PARS_LIST),256),dtype=float)
shifts = np.zeros((2,256),dtype=float) #0 = hydr, 1 = stck
 

#find initial values of the parameters from model.h file (if it gets updated, changes are read without modifying the code)
pars_from_modelh = []
vals_from_modelh = []

model_file = open("model.h",'r')
    

for line in model_file.readlines() :
    vals = line.strip().split()
    for i in range(len(PARS_LIST)):
        if len(vals) >= 3:
            if PARS_LIST[i] == vals[1] or (PARS_LIST[i]+"_OXDNA2" == vals[1] and PARS_LIST[i] != "HYDR_EPS" ) : #second condition is for FENE_R0
                pars_from_modelh.append(PARS_LIST[i])
                if vals[2]=="PI":
                    vals_from_modelh.append(math.pi)
                elif vals[2]=="(PI*0.5f)":
                    vals_from_modelh.append(math.pi*0.5)
                elif len(vals) == 5 :
                    #print(vals[2]+" " +vals[3]+" "+vals[4])
                    if vals[2]+" " +vals[3]+" "+vals[4]=="(PI - 2.35f)":
                        vals_from_modelh.append(math.pi-2.35)    
                    if vals[2]+" " +vals[3]+" "+vals[4]=="(PI - 0.875f)":
                        vals_from_modelh.append(math.pi-0.875)  
                else:
                    vals_from_modelh.append(float(vals[2][:-1]))
                break
                
        
model_file.close()

#over_indices: 
#0 - parameter index
#1 - tetramer type (convert base 4 to base 10)
#over_vals: corresponding parameter value
over_indices = []
over_vals = []
stck_fact_eps = 0.18
stck_fact_eps_read = False

#read parameters from SD file

SD_par_file = open("oxDNA_sequence_dependent_parameters_in.txt",'r')
for line in SD_par_file.readlines() :
    vals = line.strip().split()
    if len(vals) == 0:
        continue
    if vals[0] == "STCK_FACT_EPS":
        stck_fact_eps = float(vals[2])
        stck_fact_eps_read = True
        continue
    vals_cn = vals[0].split("_")
    if vals_cn[0] == "STCK" and len(vals_cn) == 3:
        for i in range(4):
            for j in range(4):
                ty=i+base_to_id(vals_cn[1])*4+base_to_id(vals_cn[2])*4*4+j*4*4*4    #from base 4 to base 10
                oi = [PARS_LIST.index("STCK_EPS"),ty]
                over_indices.append(oi)
                over_vals.append(float(vals[2]))
    elif vals_cn[0] == "HYDR" and len(vals_cn) == 3:
        for i in range(4):
            for j in range(4):
                ty=i+base_to_id(vals_cn[1])*4+base_to_id(vals_cn[2])*4*4+j*4*4*4    #from base 4 to base 10
                oi = [PARS_LIST.index("HYDR_EPS"),ty]
                over_indices.append(oi)
                over_vals.append(float(vals[2]))
    elif vals_cn[0] == "HYDR" or vals_cn[0] == "CRST" :
        par_name = vals_cn[0]
        for i in range(1,len(vals_cn)-2):
            par_name += "_"+vals_cn[i]
        for i in range(4):
            for j in range(4):
                ty=i+base_to_id(vals_cn[len(vals_cn)-2])*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4+j*4*4*4    #from base 4 to base 10
                oi = [PARS_LIST.index(par_name),ty]
                over_indices.append(oi)
                over_vals.append(float(vals[2]))
        if vals_cn[1] == "THETA2" or vals_cn[1] == "THETA5" or vals_cn[1] == "THETA7":
            par_name = vals_cn[0]
            if vals_cn[1] == "THETA2" : 
                 par_name += "_THETA3"
            elif vals_cn[1] == "THETA5" : 
                 par_name += "_THETA6"
            elif vals_cn[1] == "THETA7" : 
                 par_name += "_THETA8"
            for i in range(2,len(vals_cn)-2):
                 par_name += "_"+vals_cn[i]
                 oi = [PARS_LIST.index(par_name),ty]
                 over_indices.append(oi)
                 over_vals.append(float(vals[2]))
            for i in range(4):
                for j in range(4):
                    ty=i+base_to_id(vals_cn[len(vals_cn)-2])*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4+j*4*4*4    #from base 4 to base 10
                    oi = [PARS_LIST.index(par_name),ty]
                    over_indices.append(oi)
                    over_vals.append(float(vals[2]))




    elif vals_cn[0] == "STCK" or vals_cn[0] == "FENE" :
        par_name = vals_cn[0]
        for i in range(1,len(vals_cn)-4):
            par_name += "_"+vals_cn[i]
        ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*4+base_to_id(vals_cn[len(vals_cn)-2])*4*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4*4    #from base 4 to base 10
        oi = [PARS_LIST.index(par_name),ty]
        over_indices.append(oi)
        over_vals.append(float(vals[2]))
        if vals_cn[1] == "THETA5" :
            par_name = vals_cn[0]+"_THETA6"
            for i in range(2,len(vals_cn)-4):
                 par_name += "_"+vals_cn[i]
                 oi = [PARS_LIST.index(par_name),ty]
                 over_indices.append(oi)
                 over_vals.append(float(vals[2]))

            
SD_par_file.close()
    
if stck_fact_eps_read :
    print("STCK_FACT_EPS read from SD parameters file")
else:
    print("WARNING: No STCK_FACT_EPS found in SD parameters file")

    
#over_indices and over_vals are indices and values of parameters to overwrite
#note: overwrite pars must include STCK_x_y
#T = temperature in oxdna units. This is needed to correctly set STCK_EPS
def init_oxpars(over_indices, over_vals,T) :
    
    for i in range(len(PARS_LIST)) :
        
        if PARS_LIST[i] == "STCK_EPS":
            for j in range(256) :
                OXPS_zero[i][j] = 1.
            continue
        if PARS_LIST[i] == "HYDR_EPS":
            for j in range(256) :
                OXPS_zero[i][j] = 0.
            continue

        index = pars_from_modelh.index(PARS_LIST[i])
        val = vals_from_modelh[index]
        
        for j in range(256) :
            OXPS_zero[i][j] = val
    
    #here we use the initial custom parameters, must include stck_x_y and hydr_x_y
    if [4,48] not in over_indices :
        print("No HYDR_x_y in SD file. Terminating")
        exit(1)
    if [44,0] not in over_indices :
        print("No STCK_x_y in SD file. Terminating")
        exit(1)
    for i in range(len(over_indices)) :
        OXPS_zero[over_indices[i][0]][over_indices[i][1]] = over_vals[i]
        
    #set eps and shifts
    for j in range(256) :
        #hydr
        shifts[0][j] = Morse(OXPS_zero[par_index[7]][j],OXPS_zero[par_index[4]][j],OXPS_zero[par_index[5]][j],OXPS_zero[par_index[6]][j])
        #stacking
        OXPS_zero[par_index[44]][j] =  OXPS_zero[par_index[44]][j]* (1.0 - stck_fact_eps + (T * 9.0 * stck_fact_eps))
        shifts[1][j] = Morse(OXPS_zero[par_index[47]][j],OXPS_zero[par_index[44]][j],OXPS_zero[par_index[45]][j],OXPS_zero[par_index[46]][j])
    
    return


#print(stck_fact_eps)

T = 0.1

init_oxpars(over_indices, over_vals,T)


"""

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


PARS = copy.deepcopy(OXPS_zero) #tensor with current value of the parameterset. Use this to compute energy 
#overwrite opti parameters
for i in range(len(OPTI_P_INDICES)):
    PARS[OPTI_P_INDICES[i][0]][OPTI_P_INDICES[i][1]] = OPTI_P[i]
"""


###################################################################################################
############## READ TRAJECTORY AND COMPUTE OXDNA COORDINATES (i.e angles and distances) ###########
###################################################################################################

#define tensors

#bonded
fene_r = []
stck_r = []
th4_bn = []
th5 = []
th6 = []
cosphi1 = []
cosphi2 = []
types_bn = []

#unbonded

#compute rcut
#hydr distance cutoff. If rhydr is larger than this, then both hydr and crst are zero.

rcut_high = 0.
rcut_low = 1000.


for j in range(len(OXPS_zero[par_index[13]])):
    if OXPS_zero[par_index[13]][j] > rcut_high:
        rcut_high = OXPS_zero[par_index[13]][j]
    if OXPS_zero[par_index[85]][j] > rcut_high:
        rcut_high = OXPS_zero[par_index[85]][j]
    if OXPS_zero[par_index[124]][j] > rcut_high:
        rcut_high = OXPS_zero[par_index[124]][j]
    if OXPS_zero[par_index[12]][j] < rcut_low:
        rcut_low = OXPS_zero[par_index[12]][j]
    if OXPS_zero[par_index[84]][j] < rcut_low:
        rcut_low = OXPS_zero[par_index[84]][j]
    if OXPS_zero[par_index[123]][j] < rcut_low:
        rcut_low = OXPS_zero[par_index[123]][j]

rcut_high = rcut_high + 0.000005
rcut_low = rcut_low - 0.000005
print("rcuts: ",rcut_low,rcut_high)
rcut_sq_high = rcut_high*rcut_high
rcut_sq_low = rcut_low*rcut_low

hydr_r = []
th1 = []
th2 = []
th3 = []
th4_unbn = []
th7 = []
th8 = []
types_unbn = []


class topo :
    def __init__(self, nid, sid, bty, do, up) :
        self.id = nid
        self.strand_id = sid
        self.base_type = bty
        self.down_id = do
        self.up_id = up


pos_bb1 = [-0.34,-0.34,-0.34,-0.34]
pos_bb2 = [0.3408,0.3408,0.3408,0.3408]
#note: oxdna2
pos_stck = [0.34,0.34,0.34,0.34]
pos_hydr = [0.4,0.4,0.4,0.4]
#note: oxdna3
#pos_stck = [0.37,0.37,0.37,0.37]
#pos_hydr = [0.43,0.37,0.43,0.37]
class nucleotide :
    def __init__(self, C, BV, N, ty) :
        t = base_to_id(ty)
        self.type = t
        self.c = C
        self.bv = BV
        self.n = N
        self.norv = np.cross(N,BV)
        self.bb = C + pos_bb1[t]*BV + pos_bb2[t]*self.norv
        self.stck = C + pos_stck[t]*BV
        self.hydr = C + pos_hydr[t]*BV
    

def read_oxdna_trajectory(tr_file, topo_file):
    
    Nb = 0
    Ns = 0
    nid = 0
    
    topology = []
    counts = 0
    for line in topo_file.readlines() :
        counts+=1
        vals = line.split()
        if len(vals) == 2 :
            Nb = int(vals[0])
            Ns = int(vals[1])
            if Ns != 2 :
                print("Number of strands in not 2.")
                exit(1)
        else :
           to = topo(nid, int(vals[0]), vals[1], int(vals[2]), int(vals[3]))
           topology.append(to)
           nid += 1
           
    counts = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if a == 't':
            nid = 0
            config = []
            counts+=1
            #print(counts)
            #if counts == 501:
            #   break
            counts_b = -1
        else:
            if a != 'b' and a != 'E':
                counts_b += 1
                vals = line.split()
                c = np.array([float(vals[0]),float(vals[1]), float(vals[2])])
                bv = np.array([float(vals[3]),float(vals[4]), float(vals[5])])
                n = np.array([float(vals[6]),float(vals[7]), float(vals[8])])
                config.append(nucleotide(c,bv,n, topology[counts_b].base_type))
                nid += 1
                #configuration read. Computing coordinates
                if len(config) == Nb :
                    
                    #compute coordinates
                    fene_r_1conf = []
                    stck_r_1conf = []
                    th4_bn_1conf = []
                    th5_1conf = []
                    th6_1conf = []
                    cosphi1_1conf = []
                    cosphi2_1conf = []
                    types_1conf = []

                    #bonded basepairs
                    #find pairs and tetramer type
                    for i in range(len(topology)) :
                        ty0 = 0
                        ty1 = 0
                        ty2 = 0
                        ty3 = 0
                        if topology[i].up_id == -1: continue
                        n1 = config[i]
                        ty1 = base_to_id(topology[i].base_type)
                        n2 = n1

                        if topology[i].down_id != -1:
                         for j in range(len(topology)) :
                             if topology[j].id == topology[i].down_id:
                                ty0 = base_to_id(topology[j].base_type)
                                break
                        for j in range(len(topology)) :
                            if topology[j].id == topology[i].up_id:
                               ty2 = base_to_id(topology[j].base_type)
                               n2 = config[j]
                               if topology[j].up_id == -1:
                                   ty3 = ty0
                               else:
                                   for z in range(len(topology)) :
                                       if topology[z].id == topology[j].up_id:
                                          ty3 = base_to_id(topology[z].base_type)
                                          break
                                      
                               break
                        if topology[i].down_id == -1:
                            ty0 = ty3
                            
                        ty = ty0+ty1*4+ty2*4*4+ty3*4*4*4 #tetramer type in base 10
                        
                        types_1conf.append(ty)
                        
                        #compute bnd pair coordinates
                        
                        fene_r_1conf.append( np.linalg.norm(n1.bb-n2.bb) )
                        stck_r_1conf.append( np.linalg.norm(n1.stck-n2.stck) )
                        th4_bn_1conf.append( np.arccos(np.dot(n1.n,n2.n)) )
                        

            			#note rbb is the distance between backbone sites  in oxdna1 (see standalone)!
                        bp1 = n1.c - 0.4*n1.bv
                        bp2 = n2.c - 0.4*n2.bv
                        rbb = (bp1 - bp2)/np.linalg.norm((bp1 - bp2))
                        rstck = (n1.stck - n2.stck)/np.linalg.norm((n1.stck - n2.stck))
                        
                        th5_1conf.append( np.arccos(-np.dot(n2.n,rstck)))
                        th6_1conf.append( np.arccos(-np.dot(n1.n,rstck)))
                        cosphi1_1conf.append( np.dot(n2.norv,rbb))
                        cosphi2_1conf.append( np.dot(n1.norv,rbb))


                    types_bn.append(types_1conf)
                    fene_r.append(fene_r_1conf)
                    stck_r.append(stck_r_1conf)
                    th4_bn.append(th4_bn_1conf)
                    th5.append(th5_1conf)
                    th6.append(th6_1conf)
                    cosphi1.append(cosphi1_1conf)
                    cosphi2.append(cosphi2_1conf)
            
                    hydr_r_1conf = []
                    th1_1conf = []
                    th2_1conf = []
                    th3_1conf = []
                    th4_unbn_1conf = []
                    th7_1conf = []
                    th8_1conf = []
                    types_unbn_1conf = []
            
                    #TODO UNDBONDED
                    for i in range(len(topology)) :
                        ty0 = 0
                        ty1 = 0
                        ty2 = 0
                        ty3 = 0
                        n1 = config[i]
                        ty1 = base_to_id(topology[i].base_type)
                        for z in range(len(topology)) :
                           if topology[z].id == topology[i].down_id:
                              ty0 = base_to_id(topology[z].base_type)
                              break
                        ty0 = base_to_id(topology[i].base_type)
                        for j in range(len(topology)) :
                            n2 = config[j]
                            if np.dot(n1.hydr - n2.hydr, n1.hydr - n2.hydr) < rcut_sq_low: continue #verlet cutoff
                            if np.dot(n1.hydr - n2.hydr, n1.hydr - n2.hydr) > rcut_sq_high: continue #verlet cutoff
                            if topology[j].id <= topology[i].id: continue #ordered ids (avoid pairs repetition)
                            if topology[j].id == topology[i].down_id or topology[j].id == topology[i].up_id: continue #no bonded pairs
                            
                            ty2 = base_to_id(topology[j].base_type)
                            
                            if topology[j].up_id == -1:
                               ty3 = ty0
                            else:
                                for z in range(len(topology)) :
                                   if topology[z].id == topology[j].up_id:
                                      ty3 = base_to_id(topology[z].base_type)
                                      break
                        
                            if topology[i].down_id == -1:
                               ty0 = ty3
                               
                            ty = ty0+ty1*4+ty2*4*4+ty3*4*4*4 #tetramer type in base 10
                        
                            types_unbn_1conf.append(ty)
                        
                            #compute unbnd pair coordinates
                            hydr_r_1conf.append(np.linalg.norm((n1.hydr - n2.hydr)))
                            rhydr = (n1.hydr - n2.hydr)/np.linalg.norm((n1.hydr - n2.hydr))

                            th1_1conf.append(np.acos(-np.dot(n1.bv,n2.bv)))
                            th3_1conf.append(np.acos(-np.dot(n1.bv,rhydr)))
                            th2_1conf.append(np.acos(np.dot(n2.bv,rhydr)))
                        
                            th4_unbn_1conf.append(np.acos(np.dot(n1.n,n2.n)))
                            th8_1conf.append(np.acos(-np.dot(n1.n,rhydr)))
                            th7_1conf.append(np.acos(np.dot(n2.n,rhydr)))
                        

                    types_unbn.append(types_unbn_1conf)
                    hydr_r.append(hydr_r_1conf)

                    th1.append(th1_1conf)
                    th2.append(th2_1conf)
                    th3.append(th3_1conf)
                    
                    th4_unbn.append(th4_unbn_1conf)
                    th7.append(th7_1conf)
                    th8.append(th8_1conf)
                             
                    
    return


tr_file = open("trajectory.dat",'r')
topo_file = open("generated.top", 'r')

read_oxdna_trajectory(tr_file, topo_file)

tr_file.close()
topo_file.close()

#make unbnd tensor square. Extra unbnd pairs have zero interaction energy.

max_ints = 0
for j in range(len(types_unbn)):
    if len(types_unbn[j]) > max_ints:
       max_ints = len(types_unbn[j])
print("max unbn pairs: "+str(max_ints))
for j in range(len(types_unbn)):
    for z in range(len(types_unbn[j]),max_ints):
        types_unbn[j].append(0)
        hydr_r[j].append(0.)

        th1[j].append(0.)
        th2[j].append(0.)
        th3[j].append(0.)

        th4_unbn[j].append(0.)
        th7[j].append(0.)
        th8[j].append(0.)


#given a configuration, computes the oxdna potential energy.
#only FENE, hydrogen, stacking and cross-stacking
#takes tensors with distances and angles for every interacting pair of every configurations of multiple trajectories and computes energy with parameterset = PARS

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

FENE_R = torch.tensor(fene_r, device=device)
STCK_R = torch.tensor(stck_r, device=device)
TH4_BN = torch.tensor(th4_bn, device=device)
TH5 = torch.tensor(th5, device=device)
TH6 = torch.tensor(th6, device=device)
COSPHI1 = torch.tensor(cosphi1, device=device)
COSPHI2 = torch.tensor(cosphi2, device=device)

TYPES_BN = torch.tensor(types_bn,device=device)


HYDR_R = torch.tensor(hydr_r, device=device)
TH1 = torch.tensor(th1, device=device)
TH2 = torch.tensor(th2, device=device)
TH3 = torch.tensor(th3, device=device)
TH4_UNBN = torch.tensor(th4_unbn, device=device)
TH7 = torch.tensor(th7, device=device)
TH8 = torch.tensor(th8, device=device)

TYPES_UNBN = torch.tensor(types_unbn,device=device)

"""
pbn = np.zeros([len(PARS_LIST),len(types_bn),len(types_bn[0])])
shift_stck = np.zeros([len(types_bn),len(types_bn[0])])
shift_hydr = np.zeros([len(types_unbn),len(types_unbn[0])])
punbn = np.zeros([len(PARS_LIST),len(types_unbn),len(types_unbn[0])])


for i in range(len(PARS_LIST)):
        for j in range(len(types_bn)) :
        	for z in range(len(types_bn[j])) :
                	pbn[i][j][z] = OXPS_zero[i][types_bn[j][z]]
        for j in range(len(types_unbn)) :
        	for z in range(len(types_unbn[j])) :
                        punbn[i][j][z] = OXPS_zero[i][types_unbn[j][z]]
                    
for i in range(len(types_bn)) :
	for j in range(len(types_bn[i])) :
        	shift_stck[i][j] = shifts[1][types_bn[i][j]]

for i in range(len(types_unbn)) :
	for j in range(len(types_unbn[i])) :
        	shift_hydr[i][j] = shifts[0][types_unbn[i][j]]


PARS_BN = torch.tensor(pbn, device=device)
SHIFT_STCK = torch.tensor(shift_stck, device=device)
SHIFT_HYDR = torch.tensor(shift_hydr, device=device)
PARS_UNBN = torch.tensor(punbn, device=device)


UPDATE_MASK = torch.zeros(len(PARS_LIST),256,device=device)

UPDATE_CONSTRATINTS = torch.zeros(len(PARS_LIST),device=device)
#INDEX = torch.zeros(len(PARS_LIST))
INDEX = torch.zeros(15,len(PARS_LIST))
for i in range(len(PARS_LIST)):
    for j in range(15):
        INDEX[j][i] = i-j+8
"""

"""
PARS = torch.tensor(OXPS_zero,device=device)

for i in range(256):
    PARS[15][i] = 1.6
    PARS[20][i] = 1.4
    
PARS[25][0] = 1.6
PARS[25][2] = 1.4
#list with indices of parameters to impose continuity of
#f1
f1_r0_indices = []
f1_a_indices = []
f1_rc_indices = []
f1_bl_indices = []
f1_bh_indices = []
f1_rl_indices = []
f1_rh_indices = []
f1_rcl_indices = []
f1_rch_indices = []
#f4
f4_a_indices = []
f4_b_indices = []
f4_ts_indices = []
f4_tc_indices = []

f4_a_indices = [15,20,25]
f4_b_indices = [16,21,26]
f4_ts_indices = [17,22,27]
f4_tc_indices = [18,23,28]


f4_A = torch.tensor(f4_a_indices,device=device)
f4_B = torch.tensor(f4_b_indices,device=device)
f4_TS = torch.tensor(f4_ts_indices,device=device)
f4_TC = torch.tensor(f4_tc_indices,device=device)

PARS[f4_TS] = torch.sqrt(0.81225/PARS[f4_A])
PARS[f4_TC] = 1./PARS[f4_A]/PARS[f4_TS]
PARS[f4_B] = PARS[f4_A]*PARS[f4_TS]/(PARS[f4_TC]-PARS[f4_TS])

test_file = open("out_test_uppar.txt", 'w')
print("th1",file=test_file)
print(PARS[15],file=test_file)
print(PARS[16],file=test_file)
print(PARS[17],file=test_file)
print(PARS[18],file=test_file)

print("th2",file=test_file)
print(PARS[20],file=test_file)
print(PARS[21],file=test_file)
print(PARS[22],file=test_file)
print(PARS[23],file=test_file)

print("th3",file=test_file)
print(PARS[25],file=test_file)
print(PARS[26],file=test_file)
print(PARS[27],file=test_file)
print(PARS[28],file=test_file)

test_file.close()


#IT WORKS!!
PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 41, torch.sqrt(0.81225/PARS[INDEX[6]]),PARS[INDEX[8]])  #ts
PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 42, 1./PARS[INDEX[7]]/PARS[INDEX[5]],PARS[INDEX[8]])    #tc
PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 43, PARS[INDEX[7]]*PARS[INDEX[9]]/(PARS[INDEX[10]]-PARS[INDEX[9]]),PARS[INDEX[8]]) #b

"""
"""
def compute_energy(PARS_OPTI):
    
    #update parameters
    PARS = torch.where(UPDATE_MASK>0,PARS_OPTI,OXPS_zero)
    
    #impose symmetries
    ###todo
    
    #impose continuity
    

    #not used!
    #f2
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 21, PARS[INDEX[4]]-0.08,PARS[INDEX[8]]) #rl
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 22, PARS[INDEX[3]]-0.08,PARS[INDEX[8]]) #rh
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 23, PARS[INDEX[7]]+0.1,PARS[INDEX[8]])  #rc
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 24, PARS[INDEX[6]]-1./(PARS[INDEX[6]]-PARS[INDEX[1]])*(torch.square(PARS[INDEX[6]]-PARS[INDEX[1]])-torch.square(PARS[INDEX[2]]-PARS[INDEX[1]])),PARS[INDEX[8]])   #rcl
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 25, PARS[INDEX[6]]-1./(PARS[INDEX[6]]-PARS[INDEX[0]])*(torch.square(PARS[INDEX[6]]-PARS[INDEX[0]])-torch.square(PARS[INDEX[1]]-PARS[INDEX[0]])),PARS[INDEX[8]])   #rch
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 26, 0.5*(PARS[INDEX[6]]-PARS[INDEX[10]])/(PARS[INDEX[12]]-PARS[INDEX[10]]),PARS[INDEX[8]])  #bl
    PARS[INDEX[8]] = torch.where(UPDATE_CONSTRATINTS == 26, 0.5*(PARS[INDEX[5]]-PARS[INDEX[10]])/(PARS[INDEX[12]]-PARS[INDEX[10]]),PARS[INDEX[8]])  #bh
   
"""

""" 
for i in range(len(PARS_LIST)):
        for j in range(len(types_bn)) :
        	for z in range(len(types_bn[j])) :
                	pbn[i][j][z] = PARS[i][types_bn[j][z]]
        for j in range(len(types_unbn)) :
        	for z in range(len(types_unbn[j])) :
                        punbn[i][j][z] = PARS[i][types_unbn[j][z]]
                    
for i in range(len(types_bn)) :
	for j in range(len(types_bn[i])) :
        	shift_stck[i][j] = shifts[1][types_bn[i][j]]

for i in range(len(types_unbn)) :
	for j in range(len(types_unbn[i])) :
        	shift_hydr[i][j] = shifts[0][types_unbn[i][j]]
            
PARS_BN = torch.tensor(pbn, device=device)
SHIFT_STCK = torch.tensor(shift_stck, device=device)
SHIFT_HYDR = torch.tensor(shift_hydr, device=device)
PARS_UNBN = torch.tensor(punbn, device=device)

"""


SHIFT_STCK = torch.tensor(shifts[1], device=device)
SHIFT_HYDR = torch.tensor(shifts[0], device=device)

PARS = torch.tensor(OXPS_zero,device=device)
ZEROS_BN = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
ONES_BN = torch.ones(len(types_bn),len(types_bn[0]),device=device)
ZEROS_UNBN = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)

    
for n in range(10000) :
       
    #FENE
    energy_fene = -PARS[par_index[0]][TYPES_BN]/2.*torch.log( 1.-torch.square( FENE_R-PARS[par_index[1]][TYPES_BN] )/PARS[par_index[3]][TYPES_BN])
    
    #fene energy for each configuration                            
    energy_fene = torch.sum(energy_fene, dim=1)
    energy_fene = energy_fene/48
    
    #STACKING
    #radial part
    #piecewise energy: compute all possibilities, and filter with wher
    f1_1 = PARS[par_index[44]][TYPES_BN]*torch.square( 1.-torch.exp(-(STCK_R-PARS[par_index[45]][TYPES_BN])*PARS[par_index[46]][TYPES_BN]) )-SHIFT_STCK[TYPES_BN]
    f1_2 = PARS[par_index[44]][TYPES_BN]*PARS[par_index[48]][TYPES_BN]*torch.square( (STCK_R-PARS[par_index[52]][TYPES_BN]) )  #blow
    f1_3 = PARS[par_index[44]][TYPES_BN]*PARS[par_index[49]][TYPES_BN]*torch.square( (STCK_R-PARS[par_index[53]][TYPES_BN]) )   #bhigh
    f1_4 = ZEROS_BN
    
    #apply if conditions
    f1 = torch.where(STCK_R<PARS[par_index[53]][TYPES_BN],f1_3,f1_4)
    f1 = torch.where(STCK_R<=PARS[par_index[51]][TYPES_BN],f1_1,f1)
    f1 = torch.where(STCK_R<PARS[par_index[50]][TYPES_BN],f1_2,f1)
    f1 = torch.where(STCK_R<PARS[par_index[52]][TYPES_BN],f1_4,f1)
    
    #th4, th5, th6 Angular modulations
    
    f4_th4_1 = 1.-PARS[par_index[55]][TYPES_BN]*torch.square( (TH4_BN-PARS[par_index[54]][TYPES_BN]) )
    f4_th4_2 = PARS[par_index[56]][TYPES_BN]*torch.square( (PARS[par_index[54]][TYPES_BN]-PARS[par_index[58]][TYPES_BN]-TH4_BN) )    #low
    f4_th4_3 = PARS[par_index[56]][TYPES_BN]*torch.square( (PARS[par_index[54]][TYPES_BN]+PARS[par_index[58]][TYPES_BN]-TH4_BN) )   #high
    f4_th4_4 = ZEROS_BN
    
    f4_th4 = torch.where(TH4_BN<PARS[par_index[54]][TYPES_BN]+PARS[par_index[58]][TYPES_BN],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_BN<=PARS[par_index[54]][TYPES_BN]+PARS[par_index[57]][TYPES_BN],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_BN<PARS[par_index[54]][TYPES_BN]-PARS[par_index[57]][TYPES_BN],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_BN<PARS[par_index[54]][TYPES_BN]-PARS[par_index[58]][TYPES_BN],f4_th4_4,f4_th4)
    
    f4_th5_1 = 1.-PARS[par_index[60]][TYPES_BN]*torch.square( (TH5-PARS[par_index[59]][TYPES_BN]) )
    f4_th5_2 = PARS[par_index[61]][TYPES_BN]*torch.square( (PARS[par_index[59]][TYPES_BN]-PARS[par_index[63]][TYPES_BN]-TH5) )    #low
    f4_th5_3 = PARS[par_index[61]][TYPES_BN]*torch.square( (PARS[par_index[59]][TYPES_BN]+PARS[par_index[63]][TYPES_BN]-TH5) )   #high
    f4_th5_4 = ZEROS_BN
    
    f4_th5 = torch.where(TH5<PARS[par_index[59]][TYPES_BN]+PARS[par_index[63]][TYPES_BN],f4_th5_3,f4_th5_4)
    f4_th5 = torch.where(TH5<=PARS[par_index[59]][TYPES_BN]+PARS[par_index[62]][TYPES_BN],f4_th5_1,f4_th5)
    f4_th5 = torch.where(TH5<PARS[par_index[59]][TYPES_BN]-PARS[par_index[62]][TYPES_BN],f4_th5_2,f4_th5)
    f4_th5 = torch.where(TH5<PARS[par_index[59]][TYPES_BN]-PARS[par_index[63]][TYPES_BN],f4_th5_4,f4_th5)
    
    f4_th6_1 = 1.-PARS[par_index[65]][TYPES_BN]*torch.square( (TH6-PARS[par_index[64]][TYPES_BN]) )
    f4_th6_2 = PARS[par_index[66]][TYPES_BN]*torch.square( (PARS[par_index[64]][TYPES_BN]-PARS[par_index[68]][TYPES_BN]-TH6) )    #low
    f4_th6_3 = PARS[par_index[66]][TYPES_BN]*torch.square( (PARS[par_index[64]][TYPES_BN]+PARS[par_index[68]][TYPES_BN]-TH6) )   #high
    f4_th6_4 = ZEROS_BN
    
    f4_th6 = torch.where(TH6<PARS[par_index[64]][TYPES_BN]+PARS[par_index[68]][TYPES_BN],f4_th6_3,f4_th6_4)
    f4_th6 = torch.where(TH6<=PARS[par_index[64]][TYPES_BN]+PARS[par_index[67]][TYPES_BN],f4_th6_1,f4_th6)
    f4_th6 = torch.where(TH6<PARS[par_index[64]][TYPES_BN]-PARS[par_index[67]][TYPES_BN],f4_th6_2,f4_th6)
    f4_th6 = torch.where(TH6<PARS[par_index[64]][TYPES_BN]-PARS[par_index[68]][TYPES_BN],f4_th6_4,f4_th6)
    
    #cosphi1, cosphi2 angular modulations
    f5_phi1_1 = 1.-PARS[par_index[69]][TYPES_BN]*torch.square( COSPHI1 )
    f5_phi1_2 = PARS[par_index[70]][TYPES_BN]*torch.square( (PARS[par_index[71]]-COSPHI1 ) )    #low
    f5_phi1_3 = ONES_BN   #high
    f5_phi1_4 = ZEROS_BN
    
    
    f5_phi1 = torch.where(COSPHI1>=0,f5_phi1_3,f5_phi1_1)
    f5_phi1 = torch.where(COSPHI1<PARS[par_index[72]][TYPES_BN],f5_phi1_2,f5_phi1)
    f5_phi1 = torch.where(COSPHI1<PARS[par_index[71]][TYPES_BN],f5_phi1_4,f5_phi1)
    
    f5_phi2_1 = 1.-PARS[par_index[73]][TYPES_BN]*torch.square( COSPHI2 )
    f5_phi2_2 = PARS[par_index[74]][TYPES_BN]*torch.square( (PARS[par_index[75]][TYPES_BN]-COSPHI2 ) )    #low
    f5_phi2_3 = ONES_BN  #high
    f5_phi2_4 = ZEROS_BN
    
    
    f5_phi2 = torch.where(COSPHI2>=0,f5_phi2_3,f5_phi2_1)
    f5_phi2 = torch.where(COSPHI2<PARS[par_index[76]][TYPES_BN],f5_phi2_2,f5_phi2)
    f5_phi2 = torch.where(COSPHI2<PARS[par_index[75]][TYPES_BN],f5_phi2_4,f5_phi2)
    
    energy_stck = f1*f4_th4*f4_th5*f4_th6*f5_phi1*f5_phi2

    #stck energy for each configuration                            
    energy_stck = torch.sum(energy_stck, dim=1)
    energy_stck = energy_stck/48
     
     
    #HYDROGEN
    #radial part
    #piecewise energy: compute all possibilities, and filter with wher
    f1_1 = PARS[par_index[4]][TYPES_UNBN]*torch.square( 1.-torch.exp(-(HYDR_R-PARS[par_index[5]][TYPES_UNBN])*PARS[par_index[6]][TYPES_UNBN]) )-SHIFT_HYDR[TYPES_UNBN]
    f1_2 = PARS[par_index[4]][TYPES_UNBN]*PARS[par_index[8]][TYPES_UNBN]*torch.square( (HYDR_R-PARS[par_index[12]][TYPES_UNBN]) )  #blow
    f1_3 = PARS[par_index[4]][TYPES_UNBN]*PARS[par_index[9]][TYPES_UNBN]*torch.square( (HYDR_R-PARS[par_index[13]][TYPES_UNBN]) )   #bhigh
    f1_4 = ZEROS_UNBN
   
    #apply if conditions
    f1 = torch.where(HYDR_R<PARS[par_index[13]][TYPES_UNBN],f1_3,f1_4)
    f1 = torch.where(HYDR_R<=PARS[par_index[11]][TYPES_UNBN],f1_1,f1)
    f1 = torch.where(HYDR_R<PARS[par_index[10]][TYPES_UNBN],f1_2,f1)
    f1 = torch.where(HYDR_R<PARS[par_index[12]][TYPES_UNBN],f1_4,f1)
    
    f4_th1_1 = 1.-PARS[par_index[15]][TYPES_UNBN]*torch.square( (TH1-PARS[par_index[14]][TYPES_UNBN]) )
    f4_th1_2 = PARS[par_index[16]][TYPES_UNBN]*torch.square( (PARS[par_index[14]][TYPES_UNBN]-PARS[par_index[18]][TYPES_UNBN]-TH1) )    #low
    f4_th1_3 = PARS[par_index[16]][TYPES_UNBN]*torch.square( (PARS[par_index[14]][TYPES_UNBN]+PARS[par_index[18]][TYPES_UNBN]-TH1) )   #high
    f4_th1_4 = ZEROS_UNBN

    f4_th1 = torch.where(TH1<PARS[par_index[14]][TYPES_UNBN]+PARS[par_index[18]][TYPES_UNBN],f4_th1_3,f4_th1_4)
    f4_th1 = torch.where(TH1<=PARS[par_index[14]][TYPES_UNBN]+PARS[par_index[17]][TYPES_UNBN],f4_th1_1,f4_th1)
    f4_th1 = torch.where(TH1<PARS[par_index[14]][TYPES_UNBN]-PARS[par_index[17]][TYPES_UNBN],f4_th1_2,f4_th1)
    f4_th1 = torch.where(TH1<PARS[par_index[14]][TYPES_UNBN]-PARS[par_index[18]][TYPES_UNBN],f4_th1_4,f4_th1)
    
    f4_th2_1 = 1.-PARS[par_index[20]][TYPES_UNBN]*torch.square( (TH2-PARS[par_index[19]][TYPES_UNBN]) )
    f4_th2_2 = PARS[par_index[21]][TYPES_UNBN]*torch.square( (PARS[par_index[19]][TYPES_UNBN]-PARS[par_index[23]][TYPES_UNBN]-TH2) )    #low
    f4_th2_3 = PARS[par_index[21]][TYPES_UNBN]*torch.square( (PARS[par_index[19]][TYPES_UNBN]+PARS[par_index[23]][TYPES_UNBN]-TH2) )   #high
    f4_th2_4 = ZEROS_UNBN

    f4_th2 = torch.where(TH2<PARS[par_index[19]][TYPES_UNBN]+PARS[par_index[23]][TYPES_UNBN],f4_th2_3,f4_th2_4)
    f4_th2 = torch.where(TH2<=PARS[par_index[19]][TYPES_UNBN]+PARS[par_index[22]][TYPES_UNBN],f4_th2_1,f4_th2)
    f4_th2 = torch.where(TH2<PARS[par_index[19]][TYPES_UNBN]-PARS[par_index[22]][TYPES_UNBN],f4_th2_2,f4_th2)
    f4_th2 = torch.where(TH2<PARS[par_index[19]][TYPES_UNBN]-PARS[par_index[23]][TYPES_UNBN],f4_th2_4,f4_th2)
    
    f4_th3_1 = 1.-PARS[par_index[25]][TYPES_UNBN]*torch.square( (TH3-PARS[par_index[24]][TYPES_UNBN]) )
    f4_th3_2 = PARS[par_index[26]][TYPES_UNBN]*torch.square( (PARS[par_index[24]][TYPES_UNBN]-PARS[par_index[28]][TYPES_UNBN]-TH3) )    #low
    f4_th3_3 = PARS[par_index[26]][TYPES_UNBN]*torch.square( (PARS[par_index[24]][TYPES_UNBN]+PARS[par_index[28]][TYPES_UNBN]-TH3) )   #high
    f4_th3_4 = ZEROS_UNBN

    f4_th3 = torch.where(TH3<PARS[par_index[24]][TYPES_UNBN]+PARS[par_index[28]][TYPES_UNBN],f4_th3_3,f4_th3_4)
    f4_th3 = torch.where(TH3<=PARS[par_index[24]][TYPES_UNBN]+PARS[par_index[27]][TYPES_UNBN],f4_th3_1,f4_th3)
    f4_th3 = torch.where(TH3<PARS[par_index[24]][TYPES_UNBN]-PARS[par_index[27]][TYPES_UNBN],f4_th3_2,f4_th3)
    f4_th3 = torch.where(TH3<PARS[par_index[24]][TYPES_UNBN]-PARS[par_index[28]][TYPES_UNBN],f4_th3_4,f4_th3)
    
    f4_th4_1 = 1.-PARS[par_index[30]][TYPES_UNBN]*torch.square( (TH4_UNBN-PARS[par_index[29]][TYPES_UNBN]) )
    f4_th4_2 = PARS[par_index[31]][TYPES_UNBN]*torch.square( (PARS[par_index[29]][TYPES_UNBN]-PARS[par_index[33]][TYPES_UNBN]-TH4_UNBN) )    #low
    f4_th4_3 = PARS[par_index[31]][TYPES_UNBN]*torch.square( (PARS[par_index[29]][TYPES_UNBN]+PARS[par_index[33]][TYPES_UNBN]-TH4_UNBN) )   #high
    f4_th4_4 = ZEROS_UNBN

    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[29]][TYPES_UNBN]+PARS[par_index[33]][TYPES_UNBN],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_UNBN<=PARS[par_index[29]][TYPES_UNBN]+PARS[par_index[32]][TYPES_UNBN],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[29]][TYPES_UNBN]-PARS[par_index[32]][TYPES_UNBN],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[29]][TYPES_UNBN]-PARS[par_index[33]][TYPES_UNBN],f4_th4_4,f4_th4)
    
    f4_th7_1 = 1.-PARS[par_index[35]][TYPES_UNBN]*torch.square( (TH7-PARS[par_index[34]][TYPES_UNBN]) )
    f4_th7_2 = PARS[par_index[36]][TYPES_UNBN]*torch.square( (PARS[par_index[34]][TYPES_UNBN]-PARS[par_index[38]][TYPES_UNBN]-TH7) )    #low
    f4_th7_3 = PARS[par_index[36]][TYPES_UNBN]*torch.square( (PARS[par_index[34]][TYPES_UNBN]+PARS[par_index[38]][TYPES_UNBN]-TH7) )   #high
    f4_th7_4 = ZEROS_UNBN
    
    f4_th7 = torch.where(TH7<PARS[par_index[34]][TYPES_UNBN]+PARS[par_index[38]][TYPES_UNBN],f4_th7_3,f4_th7_4)
    f4_th7 = torch.where(TH7<=PARS[par_index[34]][TYPES_UNBN]+PARS[par_index[37]][TYPES_UNBN],f4_th7_1,f4_th7)
    f4_th7 = torch.where(TH7<PARS[par_index[34]][TYPES_UNBN]-PARS[par_index[37]][TYPES_UNBN],f4_th7_2,f4_th7)
    f4_th7 = torch.where(TH7<PARS[par_index[34]][TYPES_UNBN]-PARS[par_index[38]][TYPES_UNBN],f4_th7_4,f4_th7)
    
    f4_th8_1 = 1.-PARS[par_index[40]][TYPES_UNBN]*torch.square( (TH8-PARS[par_index[39]][TYPES_UNBN]) )
    f4_th8_2 = PARS[par_index[41]][TYPES_UNBN]*torch.square( (PARS[par_index[39]][TYPES_UNBN]-PARS[par_index[43]][TYPES_UNBN]-TH8) )    #low
    f4_th8_3 = PARS[par_index[41]][TYPES_UNBN]*torch.square( (PARS[par_index[39]][TYPES_UNBN]+PARS[par_index[43]][TYPES_UNBN]-TH8) )   #high
    f4_th8_4 = ZEROS_UNBN
    
    f4_th8 = torch.where(TH8<PARS[par_index[39]][TYPES_UNBN]+PARS[par_index[43]][TYPES_UNBN],f4_th8_3,f4_th8_4)
    f4_th8 = torch.where(TH8<=PARS[par_index[39]][TYPES_UNBN]+PARS[par_index[42]][TYPES_UNBN],f4_th8_1,f4_th8)
    f4_th8 = torch.where(TH8<PARS[par_index[39]][TYPES_UNBN]-PARS[par_index[42]][TYPES_UNBN],f4_th8_2,f4_th8)
    f4_th8 = torch.where(TH8<PARS[par_index[39]][TYPES_UNBN]-PARS[par_index[43]][TYPES_UNBN],f4_th8_4,f4_th8)
    
    energy_hydr = f1*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    
    energy_hydr = torch.sum(energy_hydr, dim=1)
    energy_hydr = energy_hydr/48

    #CROSS STACKING
    #33
    
    f2_1 = PARS[par_index[77]][TYPES_UNBN]*0.5*(torch.square( HYDR_R-PARS[par_index[78]][TYPES_UNBN])-torch.square( PARS[par_index[79]][TYPES_UNBN]-PARS[par_index[78]][TYPES_UNBN]))
    f2_2 = PARS[par_index[77]][TYPES_UNBN]*PARS[par_index[80]][TYPES_UNBN]*torch.square( (HYDR_R-PARS[par_index[84]][TYPES_UNBN]) )  #blow
    f2_3 = PARS[par_index[77]][TYPES_UNBN]*PARS[par_index[81]][TYPES_UNBN]*torch.square( (HYDR_R-PARS[par_index[85]][TYPES_UNBN]) )   #bhigh
    f2_4 = ZEROS_UNBN
   
    f2 = torch.where(HYDR_R<PARS[par_index[85]][TYPES_UNBN],f2_3,f2_4)
    f2 = torch.where(HYDR_R<=PARS[par_index[83]][TYPES_UNBN],f2_1,f2)
    f2 = torch.where(HYDR_R<PARS[par_index[82]][TYPES_UNBN],f2_2,f2)
    f2 = torch.where(HYDR_R<PARS[par_index[84]][TYPES_UNBN],f2_4,f2)
   
    f4_th1_1 = 1.-PARS[par_index[87]][TYPES_UNBN]*torch.square( (TH1-PARS[par_index[86]][TYPES_UNBN]) )
    f4_th1_2 = PARS[par_index[88]][TYPES_UNBN]*torch.square( (PARS[par_index[86]][TYPES_UNBN]-PARS[par_index[90]][TYPES_UNBN]-TH1) )    #low
    f4_th1_3 = PARS[par_index[88]][TYPES_UNBN]*torch.square( (PARS[par_index[86]][TYPES_UNBN]+PARS[par_index[90]][TYPES_UNBN]-TH1) )   #high
    f4_th1_4 = ZEROS_UNBN
   
    f4_th1 = torch.where(TH1<PARS[par_index[86]][TYPES_UNBN]+PARS[par_index[90]][TYPES_UNBN],f4_th1_3,f4_th1_4)
    f4_th1 = torch.where(TH1<=PARS[par_index[86]][TYPES_UNBN]+PARS[par_index[89]][TYPES_UNBN],f4_th1_1,f4_th1)
    f4_th1 = torch.where(TH1<PARS[par_index[86]][TYPES_UNBN]-PARS[par_index[89]][TYPES_UNBN],f4_th1_2,f4_th1)
    f4_th1 = torch.where(TH1<PARS[par_index[86]][TYPES_UNBN]-PARS[par_index[90]][TYPES_UNBN],f4_th1_4,f4_th1)
   
    f4_th2_1 = 1.-PARS[par_index[92]][TYPES_UNBN]*torch.square( (TH2-PARS[par_index[91]][TYPES_UNBN]) )
    f4_th2_2 = PARS[par_index[93]][TYPES_UNBN]*torch.square( (PARS[par_index[91]][TYPES_UNBN]-PARS[par_index[95]][TYPES_UNBN]-TH2) )    #low
    f4_th2_3 = PARS[par_index[93]][TYPES_UNBN]*torch.square( (PARS[par_index[91]][TYPES_UNBN]+PARS[par_index[95]][TYPES_UNBN]-TH2) )   #high
    f4_th2_4 = ZEROS_UNBN
   
    f4_th2 = torch.where(TH2<PARS[par_index[91]][TYPES_UNBN]+PARS[par_index[95]][TYPES_UNBN],f4_th2_3,f4_th2_4)
    f4_th2 = torch.where(TH2<=PARS[par_index[91]][TYPES_UNBN]+PARS[par_index[94]][TYPES_UNBN],f4_th2_1,f4_th2)
    f4_th2 = torch.where(TH2<PARS[par_index[91]][TYPES_UNBN]-PARS[par_index[94]][TYPES_UNBN],f4_th2_2,f4_th2)
    f4_th2 = torch.where(TH2<PARS[par_index[91]][TYPES_UNBN]-PARS[par_index[95]][TYPES_UNBN],f4_th2_4,f4_th2)
   
    f4_th3_1 = 1.-PARS[par_index[97]][TYPES_UNBN]*torch.square( (TH3-PARS[par_index[96]][TYPES_UNBN]) )
    f4_th3_2 = PARS[par_index[98]][TYPES_UNBN]*torch.square( (PARS[par_index[96]][TYPES_UNBN]-PARS[par_index[100]][TYPES_UNBN]-TH3) )    #low
    f4_th3_3 = PARS[par_index[98]][TYPES_UNBN]*torch.square( (PARS[par_index[96]][TYPES_UNBN]+PARS[par_index[100]][TYPES_UNBN]-TH3) )   #high
    f4_th3_4 = ZEROS_UNBN
   
    f4_th3 = torch.where(TH3<PARS[par_index[96]][TYPES_UNBN]+PARS[par_index[100]][TYPES_UNBN],f4_th3_3,f4_th3_4)
    f4_th3 = torch.where(TH3<=PARS[par_index[96]][TYPES_UNBN]+PARS[par_index[99]][TYPES_UNBN],f4_th3_1,f4_th3)
    f4_th3 = torch.where(TH3<PARS[par_index[96]][TYPES_UNBN]-PARS[par_index[99]][TYPES_UNBN],f4_th3_2,f4_th3)
    f4_th3 = torch.where(TH3<PARS[par_index[96]][TYPES_UNBN]-PARS[par_index[100]][TYPES_UNBN],f4_th3_4,f4_th3)
   
    f4_th4_1 = 1.-PARS[par_index[102]][TYPES_UNBN]*torch.square( (TH4_UNBN-PARS[par_index[101]][TYPES_UNBN]) )
    f4_th4_2 = PARS[par_index[103]][TYPES_UNBN]*torch.square( (PARS[par_index[101]][TYPES_UNBN]-PARS[par_index[105]][TYPES_UNBN]-TH4_UNBN) )    #low
    f4_th4_3 = PARS[par_index[103]][TYPES_UNBN]*torch.square( (PARS[par_index[101]][TYPES_UNBN]+PARS[par_index[105]][TYPES_UNBN]-TH4_UNBN) )   #high
    f4_th4_4 = ZEROS_UNBN
   
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[101]][TYPES_UNBN]+PARS[par_index[105]][TYPES_UNBN],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_UNBN<=PARS[par_index[101]][TYPES_UNBN]+PARS[par_index[104]][TYPES_UNBN],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[101]][TYPES_UNBN]-PARS[par_index[104]][TYPES_UNBN],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[101]][TYPES_UNBN]-PARS[par_index[105]][TYPES_UNBN],f4_th4_4,f4_th4)
   
    f4_th7_1 = 1.-PARS[par_index[107]][TYPES_UNBN]*torch.square( (TH7-PARS[par_index[106]][TYPES_UNBN]) )
    f4_th7_2 = PARS[par_index[108]][TYPES_UNBN]*torch.square( (PARS[par_index[106]][TYPES_UNBN]-PARS[par_index[110]][TYPES_UNBN]-TH7) )    #low
    f4_th7_3 = PARS[par_index[108]][TYPES_UNBN]*torch.square( (PARS[par_index[106]][TYPES_UNBN]+PARS[par_index[110]][TYPES_UNBN]-TH7) )   #high
    f4_th7_4 = ZEROS_UNBN
   
    f4_th7 = torch.where(TH7<PARS[par_index[106]][TYPES_UNBN]+PARS[par_index[110]][TYPES_UNBN],f4_th7_3,f4_th7_4)
    f4_th7 = torch.where(TH7<=PARS[par_index[106]][TYPES_UNBN]+PARS[par_index[109]][TYPES_UNBN],f4_th7_1,f4_th7)
    f4_th7 = torch.where(TH7<PARS[par_index[106]][TYPES_UNBN]-PARS[par_index[109]][TYPES_UNBN],f4_th7_2,f4_th7)
    f4_th7 = torch.where(TH7<PARS[par_index[106]][TYPES_UNBN]-PARS[par_index[110]][TYPES_UNBN],f4_th7_4,f4_th7)
   
    f4_th8_1 = 1.-PARS[par_index[112]][TYPES_UNBN]*torch.square( (TH8-PARS[par_index[111]][TYPES_UNBN]) )
    f4_th8_2 = PARS[par_index[113]][TYPES_UNBN]*torch.square( (PARS[par_index[111]][TYPES_UNBN]-PARS[par_index[115]][TYPES_UNBN]-TH8) )    #low
    f4_th8_3 = PARS[par_index[113]][TYPES_UNBN]*torch.square( (PARS[par_index[111]][TYPES_UNBN]+PARS[par_index[115]][TYPES_UNBN]-TH8) )   #high
    f4_th8_4 = ZEROS_UNBN
   
    f4_th8 = torch.where(TH8<PARS[par_index[111]][TYPES_UNBN]+PARS[par_index[115]][TYPES_UNBN],f4_th8_3,f4_th8_4)
    f4_th8 = torch.where(TH8<=PARS[par_index[111]][TYPES_UNBN]+PARS[par_index[114]][TYPES_UNBN],f4_th8_1,f4_th8)
    f4_th8 = torch.where(TH8<PARS[par_index[111]][TYPES_UNBN]-PARS[par_index[114]][TYPES_UNBN],f4_th8_2,f4_th8)
    f4_th8 = torch.where(TH8<PARS[par_index[111]][TYPES_UNBN]-PARS[par_index[115]][TYPES_UNBN],f4_th8_4,f4_th8)

    energy_crst = f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    
    #55
    f2_1 = PARS[par_index[116]][TYPES_UNBN]*0.5*(torch.square( HYDR_R-PARS[par_index[117]][TYPES_UNBN])-torch.square( PARS[par_index[118]][TYPES_UNBN]-PARS[par_index[117]][TYPES_UNBN]))
    f2_2 = PARS[par_index[116]][TYPES_UNBN]*PARS[par_index[119]][TYPES_UNBN]*torch.square( (HYDR_R-PARS[par_index[123]][TYPES_UNBN]) )  #blow
    f2_3 = PARS[par_index[116]][TYPES_UNBN]*PARS[par_index[120]][TYPES_UNBN]*torch.square( (HYDR_R-PARS[par_index[124]][TYPES_UNBN]) )   #bhigh
    f2_4 = ZEROS_UNBN
    
    f2 = torch.where(HYDR_R<PARS[par_index[124]][TYPES_UNBN],f2_3,f2_4)
    f2 = torch.where(HYDR_R<=PARS[par_index[122]][TYPES_UNBN],f2_1,f2)
    f2 = torch.where(HYDR_R<PARS[par_index[121]][TYPES_UNBN],f2_2,f2)
    f2 = torch.where(HYDR_R<PARS[par_index[123]][TYPES_UNBN],f2_4,f2)
    
    f4_th1_1 = 1.-PARS[par_index[126]][TYPES_UNBN]*torch.square( (TH1-PARS[par_index[125]][TYPES_UNBN]) )
    f4_th1_2 = PARS[par_index[127]][TYPES_UNBN]*torch.square( (PARS[par_index[125]][TYPES_UNBN]-PARS[par_index[129]][TYPES_UNBN]-TH1) )    #low
    f4_th1_3 = PARS[par_index[127]][TYPES_UNBN]*torch.square( (PARS[par_index[125]][TYPES_UNBN]+PARS[par_index[129]][TYPES_UNBN]-TH1) )   #high
    f4_th1_4 = ZEROS_UNBN
    
    f4_th1 = torch.where(TH1<PARS[par_index[125]][TYPES_UNBN]+PARS[par_index[129]][TYPES_UNBN],f4_th1_3,f4_th1_4)
    f4_th1 = torch.where(TH1<=PARS[par_index[125]][TYPES_UNBN]+PARS[par_index[128]][TYPES_UNBN],f4_th1_1,f4_th1)
    f4_th1 = torch.where(TH1<PARS[par_index[125]][TYPES_UNBN]-PARS[par_index[128]][TYPES_UNBN],f4_th1_2,f4_th1)
    f4_th1 = torch.where(TH1<PARS[par_index[125]][TYPES_UNBN]-PARS[par_index[129]][TYPES_UNBN],f4_th1_4,f4_th1)
    
    f4_th2_1 = 1.-PARS[par_index[131]][TYPES_UNBN]*torch.square( (TH2-PARS[par_index[130]][TYPES_UNBN]) )
    f4_th2_2 = PARS[par_index[132]][TYPES_UNBN]*torch.square( (PARS[par_index[130]][TYPES_UNBN]-PARS[par_index[134]][TYPES_UNBN]-TH2) )    #low
    f4_th2_3 = PARS[par_index[132]][TYPES_UNBN]*torch.square( (PARS[par_index[130]][TYPES_UNBN]+PARS[par_index[134]][TYPES_UNBN]-TH2) )   #high
    f4_th2_4 = ZEROS_UNBN
    
    f4_th2 = torch.where(TH2<PARS[par_index[130]][TYPES_UNBN]+PARS[par_index[134]][TYPES_UNBN],f4_th2_3,f4_th2_4)
    f4_th2 = torch.where(TH2<=PARS[par_index[130]][TYPES_UNBN]+PARS[par_index[133]][TYPES_UNBN],f4_th2_1,f4_th2)
    f4_th2 = torch.where(TH2<PARS[par_index[130]][TYPES_UNBN]-PARS[par_index[133]][TYPES_UNBN],f4_th2_2,f4_th2)
    f4_th2 = torch.where(TH2<PARS[par_index[130]][TYPES_UNBN]-PARS[par_index[134]][TYPES_UNBN],f4_th2_4,f4_th2)
    
    f4_th3_1 = 1.-PARS[par_index[136]][TYPES_UNBN]*torch.square( (TH3-PARS[par_index[135]][TYPES_UNBN]) )
    f4_th3_2 = PARS[par_index[137]][TYPES_UNBN]*torch.square( (PARS[par_index[135]][TYPES_UNBN]-PARS[par_index[139]][TYPES_UNBN]-TH3) )    #low
    f4_th3_3 = PARS[par_index[137]][TYPES_UNBN]*torch.square( (PARS[par_index[135]][TYPES_UNBN]+PARS[par_index[139]][TYPES_UNBN]-TH3) )   #high
    f4_th3_4 = ZEROS_UNBN
    
    f4_th3 = torch.where(TH3<PARS[par_index[135]][TYPES_UNBN]+PARS[par_index[139]][TYPES_UNBN],f4_th3_3,f4_th3_4)
    f4_th3 = torch.where(TH3<=PARS[par_index[135]][TYPES_UNBN]+PARS[par_index[138]][TYPES_UNBN],f4_th3_1,f4_th3)
    f4_th3 = torch.where(TH3<PARS[par_index[135]][TYPES_UNBN]-PARS[par_index[138]][TYPES_UNBN],f4_th3_2,f4_th3)
    f4_th3 = torch.where(TH3<PARS[par_index[135]][TYPES_UNBN]-PARS[par_index[139]][TYPES_UNBN],f4_th3_4,f4_th3)
    
    f4_th4_1 = 1.-PARS[par_index[141]][TYPES_UNBN]*torch.square( (TH4_UNBN-PARS[par_index[140]][TYPES_UNBN]) )
    f4_th4_2 = PARS[par_index[142]][TYPES_UNBN]*torch.square( (PARS[par_index[140]][TYPES_UNBN]-PARS[par_index[144]][TYPES_UNBN]-TH4_UNBN) )    #low
    f4_th4_3 = PARS[par_index[142]][TYPES_UNBN]*torch.square( (PARS[par_index[140]][TYPES_UNBN]+PARS[par_index[144]][TYPES_UNBN]-TH4_UNBN) )   #high
    f4_th4_4 = ZEROS_UNBN
    
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[140]][TYPES_UNBN]+PARS[par_index[144]][TYPES_UNBN],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_UNBN<=PARS[par_index[140]][TYPES_UNBN]+PARS[par_index[143]][TYPES_UNBN],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[140]][TYPES_UNBN]-PARS[par_index[143]][TYPES_UNBN],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS[par_index[140]][TYPES_UNBN]-PARS[par_index[144]][TYPES_UNBN],f4_th4_4,f4_th4)
    
    f4_th7_1 = 1.-PARS[par_index[146]][TYPES_UNBN]*torch.square( (TH7-PARS[par_index[145]][TYPES_UNBN]) )
    f4_th7_2 = PARS[par_index[147]][TYPES_UNBN]*torch.square( (PARS[par_index[145]][TYPES_UNBN]-PARS[par_index[149]][TYPES_UNBN]-TH7) )    #low
    f4_th7_3 = PARS[par_index[147]][TYPES_UNBN]*torch.square( (PARS[par_index[145]][TYPES_UNBN]+PARS[par_index[149]][TYPES_UNBN]-TH7) )   #high
    f4_th7_4 = ZEROS_UNBN
    
    f4_th7 = torch.where(TH7<PARS[par_index[145]][TYPES_UNBN]+PARS[par_index[149]][TYPES_UNBN],f4_th7_3,f4_th7_4)
    f4_th7 = torch.where(TH7<=PARS[par_index[145]][TYPES_UNBN]+PARS[par_index[148]][TYPES_UNBN],f4_th7_1,f4_th7)
    f4_th7 = torch.where(TH7<PARS[par_index[145]][TYPES_UNBN]-PARS[par_index[148]][TYPES_UNBN],f4_th7_2,f4_th7)
    f4_th7 = torch.where(TH7<PARS[par_index[145]][TYPES_UNBN]-PARS[par_index[149]][TYPES_UNBN],f4_th7_4,f4_th7)
    
    f4_th8_1 = 1.-PARS[par_index[151]][TYPES_UNBN]*torch.square( (TH8-PARS[par_index[150]][TYPES_UNBN]) )
    f4_th8_2 = PARS[par_index[152]][TYPES_UNBN]*torch.square( (PARS[par_index[150]][TYPES_UNBN]-PARS[par_index[154]][TYPES_UNBN]-TH8) )    #low
    f4_th8_3 = PARS[par_index[152]][TYPES_UNBN]*torch.square( (PARS[par_index[150]][TYPES_UNBN]+PARS[par_index[154]][TYPES_UNBN]-TH8) )   #high
    f4_th8_4 = ZEROS_UNBN

    f4_th8 = torch.where(TH8<PARS[par_index[150]][TYPES_UNBN]+PARS[par_index[154]][TYPES_UNBN],f4_th8_3,f4_th8_4)
    f4_th8 = torch.where(TH8<=PARS[par_index[150]][TYPES_UNBN]+PARS[par_index[153]][TYPES_UNBN],f4_th8_1,f4_th8)
    f4_th8 = torch.where(TH8<PARS[par_index[150]][TYPES_UNBN]-PARS[par_index[153]][TYPES_UNBN],f4_th8_2,f4_th8)
    f4_th8 = torch.where(TH8<PARS[par_index[150]][TYPES_UNBN]-PARS[par_index[154]][TYPES_UNBN],f4_th8_4,f4_th8)
    
    energy_crst = energy_crst +  f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8

    energy_crst = torch.sum(energy_crst, dim=1)
    energy_crst = energy_crst/48

    energy_tot = energy_fene + energy_stck + energy_hydr + energy_crst
#print(counts)

print("Energy_fene in:", energy_fene.device)

test_file = open("out_test.txt", 'w')
print(energy_fene, file=test_file)

print(energy_stck, file=test_file)

print(energy_hydr, file=test_file)

print(energy_crst, file=test_file)

print(energy_tot, file=test_file)

test_file.close()


