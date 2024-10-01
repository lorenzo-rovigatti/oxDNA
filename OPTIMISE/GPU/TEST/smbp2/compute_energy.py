
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
"""
out_file = open("out_test.txt",'w')
for i in range(len(pars_from_modelh)) :
    print(pars_from_modelh[i]+" "+ str(vals_from_modelh[i]),file=out_file)
out_file.close()
"""

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
    
"""
out_test1 = open("out_test1.txt",'w')

for i in range(len(over_indices)) :
    ty0 = over_indices[i][1]%4
    ty1 = (over_indices[i][1]//4)%4
    ty2 = (over_indices[i][1]//4//4)%4
    ty3 = (over_indices[i][1]//4//4//4)%4
    print(str(over_indices[i][0])+" "+str(over_indices[i][1])+" "+PARS_LIST[over_indices[i][0]]+"_"+bases[ty0]+"_"+bases[ty1]+"_"+bases[ty2]+"_"+bases[ty3] + " = " + str(over_vals[i]), file = out_test1)
    
out_test1.close()
"""

    
#over_indices and over_vals are indices and values of parameters to overwrite
#note: overwrite pars must include STCK_x_y
#T = temperature in oxdna units. This is needed to correctly set STCK_EPS
def init_oxpars(over_indices, over_vals,T) :
    
    for i in range(len(PARS_LIST)) :
        
        if PARS_LIST[i] == "STCK_EPS":
            for j in range(256) :
                OXPS_zero[i][j] = 1.
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
out_test = open("out_test2.txt",'w')
print(OXPS_zero, file=out_test)
out_test.close()
"""

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
#norvs = []

#unbonded

#compute rcut
#hydr distance cutoff. If rhydr is larger than this, then both hydr and crst are zero.

rcut_high = 0.
rcut_low = 1000.

for j in range(len(OXPS_zero[par_index[13]])):
    if OXPS_zero[par_index[13]][j] > rcut_high:
        rcut = OXPS_zero[par_index[13]][j]
    if OXPS_zero[par_index[85]][j] > rcut_high:
        rcut = OXPS_zero[par_index[85]][j]
    if OXPS_zero[par_index[124]][j] > rcut_high:
        rcut = OXPS_zero[par_index[124]][j]
    if OXPS_zero[par_index[12]][j] < rcut_low:
        rcut = OXPS_zero[par_index[12]][j]
    if OXPS_zero[par_index[84]][j] < rcut_low:
        rcut = OXPS_zero[par_index[84]][j]
    if OXPS_zero[par_index[123]][j] < rcut_low:
        rcut = OXPS_zero[par_index[123]][j]

rcut_high = rcut_high + 0.000005
rcut_low = rcut_low - 0.000005
print("rcuts: ",rcut_low,rcut_high)
rcut_sq_high = rcut_high*rcut_high
rcut_sq_low = rcut_high*rcut_low

hydr_r = []
th1 = []
th2 = []
th3 = []
th4_unbn = []
th7 = []
th8 = []
types_unbn = []
ids_list = []


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
    
    #for i in range(len(topology)) :
    #    print(topology[i].id, topology[i].base_type, topology[i].down_id,topology[i].up_id)
           
    counts = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if a == 't':
            nid = 0
            #print(a)
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
                #print(vals[0])
                c = np.array([float(vals[0]),float(vals[1]), float(vals[2])])
                #print(c)
                bv = np.array([float(vals[3]),float(vals[4]), float(vals[5])])
                #print(bv)
                n = np.array([float(vals[6]),float(vals[7]), float(vals[8])])
                #print(n)
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
                    #norv_1conf = []
                    
                    #print(counts,counts_b)
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
                        
                        #print(bases[ty0], bases[ty1], bases[ty2], bases[ty3])
                        
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
                        #norv_1conf.append(n1.norv)
                        #norv_1conf.append(rbb*np.linalg.norm((bp1 - bp2)))

                    types_bn.append(types_1conf)
                    fene_r.append(fene_r_1conf)
                    stck_r.append(stck_r_1conf)
                    th4_bn.append(th4_bn_1conf)
                    th5.append(th5_1conf)
                    th6.append(th6_1conf)
                    cosphi1.append(cosphi1_1conf)
                    cosphi2.append(cosphi2_1conf)
                    #norvs.append(norv_1conf)
            
                    hydr_r_1conf = []
                    th1_1conf = []
                    th2_1conf = []
                    th3_1conf = []
                    th4_unbn_1conf = []
                    th7_1conf = []
                    th8_1conf = []
                    types_unbn_1conf = []
                    
                    ids_1conf = []
            
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
                               
                               
                            IDS = [topology[i].id,topology[j].id]
                            ids_1conf.append(IDS)
                               
                            ty = ty0+ty1*4+ty2*4*4+ty3*4*4*4 #tetramer type in base 10
                        
                            #print(bases[ty0], bases[ty1], bases[ty2], bases[ty3])
                        
                            types_unbn_1conf.append(ty)
                        
                            #compute unbnd pair coordinates
                            hydr_r_1conf.append(np.linalg.norm((n1.hydr - n2.hydr)))
                            rhydr = (n1.hydr - n2.hydr)/np.linalg.norm((n1.hydr - n2.hydr))

                            th1_1conf.append(np.acos(-np.dot(n1.bv,n2.bv)))
                            th2_1conf.append(np.acos(-np.dot(n1.bv,rhydr)))
                            th3_1conf.append(np.acos(np.dot(n2.bv,rhydr)))
                        
                            th4_unbn_1conf.append(np.acos(-np.dot(n1.n,n2.n)))
                            th7_1conf.append(np.acos(-np.dot(n1.n,rhydr)))
                            th8_1conf.append(np.acos(np.dot(n2.n,rhydr)))
                        

                    types_unbn.append(types_unbn_1conf)
                    hydr_r.append(hydr_r_1conf)

                    th1.append(th1_1conf)
                    th2.append(th2_1conf)
                    th3.append(th3_1conf)
                    
                    th4_unbn.append(th4_unbn_1conf)
                    th7.append(th7_1conf)
                    th8.append(th8_1conf)
                    
                    ids_list.append(ids_1conf)
                        
                             
                    
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

        ids_list[j].append([-1,-1])

"""
test_file = open("out_test3.txt",'w')
print(len(FENE_R))
print(len(FENE_R[0]))
print(FENE_R, file=test_file)
print(TYPES_BN, file=test_file)
test_file.close()
"""

#given a configuration, computes the oxdna potential energy.
#only FENE, hydrogen, stacking and cross-stacking
#takes tensors with distances and angles for every interacting pair of every configurations of multiple trajectories and computes energy with parameterset = PARS



#generate tensors with parameters. 1 for bonded pairs, 1 for unbonded pairs
pbn = np.zeros([len(PARS_LIST),len(types_bn),len(types_bn[0])])
shift_stck = np.zeros([len(types_bn),len(types_bn[0])])
shift_hydr = np.zeros([len(types_unbn),len(types_unbn[0])])
punbn = np.zeros([len(PARS_LIST),len(types_unbn),len(types_unbn[0])])


#print(OXPS_zero[0])
#print(TYPES_BN[0])
for i in range(len(PARS_LIST)):
        for j in range(len(types_bn)) :
        	for z in range(len(types_bn[j])) :
                	pbn[i][j][z] = OXPS_zero[i][types_bn[j][z]]
        for j in range(len(types_unbn)) :
        	for z in range(len(types_unbn[j])) :
                        punbn[i][j][z] = OXPS_zero[i][types_unbn[j][z]]
                    
for i in range(len(types_bn)) :
	for j in range(len(types_bn[j])) :
        	shift_stck[i][j] = shifts[1][types_bn[i][j]]

for i in range(len(types_unbn)) :
	for j in range(len(types_unbn[j])) :
        	shift_hydr[i][j] = shifts[0][types_unbn[i][j]]


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

PARS_BN = torch.tensor(pbn, device=device)
SHIFT_STCK = torch.tensor(shift_stck, device=device)
SHIFT_HYDR = torch.tensor(shift_hydr, device=device)
PARS_UNBN = torch.tensor(punbn, device=device)

FENE_R = torch.tensor(fene_r, device=device)
STCK_R = torch.tensor(stck_r, device=device)
TH4_BN = torch.tensor(th4_bn, device=device)
TH5 = torch.tensor(th5, device=device)
TH6 = torch.tensor(th6, device=device)
COSPHI1 = torch.tensor(cosphi1, device=device)
COSPHI2 = torch.tensor(cosphi2, device=device)
#NORVS = torch.tensor(norvs, device=device)
#fene_r = torch.tensor(FENE_R)
#PARS_BN = torch.tensor(pbn)

HYDR_R = torch.tensor(hydr_r, device=device)
TH1 = torch.tensor(th1, device=device)
TH2 = torch.tensor(th2, device=device)
TH3 = torch.tensor(th3, device=device)
TH4_UNBN = torch.tensor(th4_unbn, device=device)
TH7 = torch.tensor(th7, device=device)
TH8 = torch.tensor(th8, device=device)
IDS_LIST = torch.tensor(ids_list, device=device)


#print("fene_r in:", fene_r.device)
#print("PARS_BN in:", PARS_BN.device)

test_file = open("out_test1.txt",'w')
#print(PARS_BN[par_index[64]],file=test_file)
#print(PARS_BN[par_index[65]],file=test_file)
#print(PARS_BN[par_index[66]],file=test_file)
#print(PARS_BN[par_index[67]],file=test_file)
#print(PARS_BN[par_index[44]],file=test_file)
#print(SHIFT_STCK,file=test_file)
#print(STCK_R,file=test_file)
#print(TH4_BN,file=test_file)
#print(TH5,file=test_file)
#print(TH6,file=test_file)
#print(COSPHI1,file=test_file)
#print(COSPHI2,file=test_file)
#print(NORVS,file=test_file)
print(IDS_LIST,file=test_file)
print(HYDR_R,file=test_file)
print(TH1,file=test_file)
print(TH2,file=test_file)
print(TH3,file=test_file)
print(TH4_UNBN,file=test_file)
print(TH7,file=test_file)
print(TH8,file=test_file)
test_file.close()


#oxp = torch.tensor(OXPS_zero)

#fene energy for each pair

for n in range(1) :

    #FENE
    energy_fene = -PARS_BN[par_index[0]]/2.*torch.log( 1.-torch.square( FENE_R-PARS_BN[par_index[1]] )/PARS_BN[par_index[3]])
    
    #fene energy for each configuration                            
    energy_fene = torch.sum(energy_fene, dim=1)
    energy_fene = energy_fene/48
    
    #STACKING
    #radial part
    #piecewise energy: compute all possibilities, and filter with wher
    f1_1 = PARS_BN[par_index[44]]*torch.square( 1.-torch.exp(-(STCK_R-PARS_BN[par_index[45]])*PARS_BN[par_index[46]]) )-SHIFT_STCK
    f1_2 = PARS_BN[par_index[44]]*PARS_BN[par_index[48]]*torch.square( (STCK_R-PARS_BN[par_index[52]]) )  #blow
    f1_3 = PARS_BN[par_index[44]]*PARS_BN[par_index[49]]*torch.square( (STCK_R-PARS_BN[par_index[53]]) )   #bhigh
    f1_4 = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
    
    #apply if conditions
    f1 = torch.where(STCK_R<PARS_BN[par_index[53]],f1_3,f1_4)
    f1 = torch.where(STCK_R<=PARS_BN[par_index[51]],f1_1,f1)
    f1 = torch.where(STCK_R<PARS_BN[par_index[50]],f1_2,f1)
    f1 = torch.where(STCK_R<PARS_BN[par_index[52]],f1_4,f1)
    
    #th4, th5, th6 Angular modulations
    
    f4_th4_1 = 1.-PARS_BN[par_index[55]]*torch.square( (TH4_BN-PARS_BN[par_index[54]]) )
    f4_th4_2 = PARS_BN[par_index[56]]*torch.square( (PARS_BN[par_index[54]]-PARS_BN[par_index[58]]-TH4_BN) )    #low
    f4_th4_3 = PARS_BN[par_index[56]]*torch.square( (PARS_BN[par_index[54]]+PARS_BN[par_index[58]]-TH4_BN) )   #high
    f4_th4_4 = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
    
    f4_th4 = torch.where(TH4_BN<PARS_BN[par_index[54]]+PARS_BN[par_index[58]],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_BN<=PARS_BN[par_index[54]]+PARS_BN[par_index[57]],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_BN<PARS_BN[par_index[54]]-PARS_BN[par_index[57]],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_BN<PARS_BN[par_index[54]]-PARS_BN[par_index[58]],f4_th4_4,f4_th4)
    
    f4_th5_1 = 1.-PARS_BN[par_index[60]]*torch.square( (TH5-PARS_BN[par_index[59]]) )
    f4_th5_2 = PARS_BN[par_index[61]]*torch.square( (PARS_BN[par_index[59]]-PARS_BN[par_index[63]]-TH5) )    #low
    f4_th5_3 = PARS_BN[par_index[61]]*torch.square( (PARS_BN[par_index[59]]+PARS_BN[par_index[63]]-TH5) )   #high
    f4_th5_4 = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
    
    f4_th5 = torch.where(TH5<PARS_BN[par_index[59]]+PARS_BN[par_index[63]],f4_th5_3,f4_th5_4)
    f4_th5 = torch.where(TH5<=PARS_BN[par_index[59]]+PARS_BN[par_index[62]],f4_th5_1,f4_th5)
    f4_th5 = torch.where(TH5<PARS_BN[par_index[59]]-PARS_BN[par_index[62]],f4_th5_2,f4_th5)
    f4_th5 = torch.where(TH5<PARS_BN[par_index[59]]-PARS_BN[par_index[63]],f4_th5_4,f4_th5)
    
    f4_th6_1 = 1.-PARS_BN[par_index[65]]*torch.square( (TH6-PARS_BN[par_index[64]]) )
    f4_th6_2 = PARS_BN[par_index[66]]*torch.square( (PARS_BN[par_index[64]]-PARS_BN[par_index[68]]-TH6) )    #low
    f4_th6_3 = PARS_BN[par_index[66]]*torch.square( (PARS_BN[par_index[64]]+PARS_BN[par_index[68]]-TH6) )   #high
    f4_th6_4 = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
    
    f4_th6 = torch.where(TH6<PARS_BN[par_index[64]]+PARS_BN[par_index[68]],f4_th6_3,f4_th6_4)
    f4_th6 = torch.where(TH6<=PARS_BN[par_index[64]]+PARS_BN[par_index[67]],f4_th6_1,f4_th6)
    f4_th6 = torch.where(TH6<PARS_BN[par_index[64]]-PARS_BN[par_index[67]],f4_th6_2,f4_th6)
    f4_th6 = torch.where(TH6<PARS_BN[par_index[64]]-PARS_BN[par_index[68]],f4_th6_4,f4_th6)
    
    #cosphi1, cosphi2 angular modulations
    f5_phi1_1 = 1.-PARS_BN[par_index[69]]*torch.square( COSPHI1 )
    f5_phi1_2 = PARS_BN[par_index[70]]*torch.square( (PARS_BN[par_index[71]]-COSPHI1 ) )    #low
    f5_phi1_3 = torch.ones(len(types_bn),len(types_bn[0]),device=device)   #high
    f5_phi1_4 = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
    
    
    f5_phi1 = torch.where(COSPHI1>=0,f5_phi1_3,f5_phi1_1)
    f5_phi1 = torch.where(COSPHI1<PARS_BN[par_index[72]],f5_phi1_2,f5_phi1)
    f5_phi1 = torch.where(COSPHI1<PARS_BN[par_index[71]],f5_phi1_4,f5_phi1)
    
    f5_phi2_1 = 1.-PARS_BN[par_index[73]]*torch.square( COSPHI2 )
    f5_phi2_2 = PARS_BN[par_index[74]]*torch.square( (PARS_BN[par_index[75]]-COSPHI2 ) )    #low
    f5_phi2_3 = torch.ones(len(types_bn),len(types_bn[0]),device=device)   #high
    f5_phi2_4 = torch.zeros(len(types_bn),len(types_bn[0]),device=device)
    
    
    f5_phi2 = torch.where(COSPHI2>=0,f5_phi2_3,f5_phi2_1)
    f5_phi2 = torch.where(COSPHI2<PARS_BN[par_index[76]],f5_phi2_2,f5_phi2)
    f5_phi2 = torch.where(COSPHI2<PARS_BN[par_index[75]],f5_phi2_4,f5_phi2)
    
    energy_stck = f1*f4_th4*f4_th5*f4_th6*f5_phi1*f5_phi2

    #stck energy for each configuration                            
    energy_stck = torch.sum(energy_stck, dim=1)
    energy_stck = energy_stck/48
    
""" 
print("Energy_fene in:", energy_fene.device)

test_file = open("out_test.txt", 'w')
print(energy_fene, file=test_file)

print("Energy_fene in:", energy_stck.device)

print(energy_stck, file=test_file)

print(f1, file=test_file)
print(f4_th4, file=test_file)
print(f4_th5, file=test_file)
print(f4_th6, file=test_file)
print(f5_phi1, file=test_file)
print(f5_phi2, file=test_file)

test_file.close()
"""
    
    
for n in range(1) :    
    #HYDROGEN
    #radial part
    #piecewise energy: compute all possibilities, and filter with wher
    f1_1 = PARS_UNBN[par_index[4]]*torch.square( 1.-torch.exp(-(HYDR_R-PARS_UNBN[par_index[5]])*PARS_UNBN[par_index[6]]) )-SHIFT_HYDR
    f1_2 = PARS_UNBN[par_index[4]]*PARS_UNBN[par_index[8]]*torch.square( (HYDR_R-PARS_UNBN[par_index[12]]) )  #blow
    f1_3 = PARS_UNBN[par_index[4]]*PARS_UNBN[par_index[9]]*torch.square( (HYDR_R-PARS_UNBN[par_index[13]]) )   #bhigh
    f1_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    #apply if conditions
    f1 = torch.where(HYDR_R<PARS_UNBN[par_index[13]],f1_3,f1_4)
    f1 = torch.where(HYDR_R<=PARS_UNBN[par_index[11]],f1_1,f1)
    f1 = torch.where(HYDR_R<PARS_UNBN[par_index[10]],f1_2,f1)
    f1 = torch.where(HYDR_R<PARS_UNBN[par_index[12]],f1_4,f1)
    
    f4_th1_1 = 1.-PARS_UNBN[par_index[15]]*torch.square( (TH1-PARS_UNBN[par_index[14]]) )
    f4_th1_2 = PARS_UNBN[par_index[16]]*torch.square( (PARS_UNBN[par_index[14]]-PARS_UNBN[par_index[18]]-TH1) )    #low
    f4_th1_3 = PARS_UNBN[par_index[16]]*torch.square( (PARS_UNBN[par_index[14]]+PARS_UNBN[par_index[18]]-TH1) )   #high
    f4_th1_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)

    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[14]]+PARS_UNBN[par_index[18]],f4_th1_3,f4_th1_4)
    f4_th1 = torch.where(TH1<=PARS_UNBN[par_index[14]]+PARS_UNBN[par_index[17]],f4_th1_1,f4_th1)
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[14]]-PARS_UNBN[par_index[17]],f4_th1_2,f4_th1)
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[14]]-PARS_UNBN[par_index[18]],f4_th1_4,f4_th1)
    
    f4_th2_1 = 1.-PARS_UNBN[par_index[20]]*torch.square( (TH2-PARS_UNBN[par_index[19]]) )
    f4_th2_2 = PARS_UNBN[par_index[21]]*torch.square( (PARS_UNBN[par_index[19]]-PARS_UNBN[par_index[23]]-TH2) )    #low
    f4_th2_3 = PARS_UNBN[par_index[21]]*torch.square( (PARS_UNBN[par_index[19]]+PARS_UNBN[par_index[23]]-TH2) )   #high
    f4_th2_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)

    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[19]]+PARS_UNBN[par_index[23]],f4_th2_3,f4_th2_4)
    f4_th2 = torch.where(TH2<=PARS_UNBN[par_index[19]]+PARS_UNBN[par_index[22]],f4_th2_1,f4_th2)
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[19]]-PARS_UNBN[par_index[22]],f4_th2_2,f4_th2)
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[19]]-PARS_UNBN[par_index[23]],f4_th2_4,f4_th2)
    
    f4_th3_1 = 1.-PARS_UNBN[par_index[25]]*torch.square( (TH3-PARS_UNBN[par_index[24]]) )
    f4_th3_2 = PARS_UNBN[par_index[26]]*torch.square( (PARS_UNBN[par_index[24]]-PARS_UNBN[par_index[28]]-TH3) )    #low
    f4_th3_3 = PARS_UNBN[par_index[26]]*torch.square( (PARS_UNBN[par_index[24]]+PARS_UNBN[par_index[28]]-TH3) )   #high
    f4_th3_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)

    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[24]]+PARS_UNBN[par_index[28]],f4_th3_3,f4_th3_4)
    f4_th3 = torch.where(TH3<=PARS_UNBN[par_index[24]]+PARS_UNBN[par_index[27]],f4_th3_1,f4_th3)
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[24]]-PARS_UNBN[par_index[27]],f4_th3_2,f4_th3)
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[24]]-PARS_UNBN[par_index[28]],f4_th3_4,f4_th3)
    
    f4_th4_1 = 1.-PARS_UNBN[par_index[30]]*torch.square( (TH4_UNBN-PARS_UNBN[par_index[29]]) )
    f4_th4_2 = PARS_UNBN[par_index[31]]*torch.square( (PARS_UNBN[par_index[29]]-PARS_UNBN[par_index[33]]-TH4_UNBN) )    #low
    f4_th4_3 = PARS_UNBN[par_index[31]]*torch.square( (PARS_UNBN[par_index[29]]+PARS_UNBN[par_index[33]]-TH4_UNBN) )   #high
    f4_th4_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)

    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[29]]+PARS_UNBN[par_index[33]],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_UNBN<=PARS_UNBN[par_index[29]]+PARS_UNBN[par_index[32]],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[29]]-PARS_UNBN[par_index[32]],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[29]]-PARS_UNBN[par_index[33]],f4_th4_4,f4_th4)
    
    f4_th7_1 = 1.-PARS_UNBN[par_index[35]]*torch.square( (TH7-PARS_UNBN[par_index[34]]) )
    f4_th7_2 = PARS_UNBN[par_index[36]]*torch.square( (PARS_UNBN[par_index[34]]-PARS_UNBN[par_index[38]]-TH7) )    #low
    f4_th7_3 = PARS_UNBN[par_index[36]]*torch.square( (PARS_UNBN[par_index[34]]+PARS_UNBN[par_index[38]]-TH7) )   #high
    f4_th7_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[34]]+PARS_UNBN[par_index[38]],f4_th7_3,f4_th7_4)
    f4_th7 = torch.where(TH7<=PARS_UNBN[par_index[34]]+PARS_UNBN[par_index[37]],f4_th7_1,f4_th7)
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[34]]-PARS_UNBN[par_index[37]],f4_th7_2,f4_th7)
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[34]]-PARS_UNBN[par_index[38]],f4_th7_4,f4_th7)
    
    f4_th8_1 = 1.-PARS_UNBN[par_index[40]]*torch.square( (TH8-PARS_UNBN[par_index[39]]) )
    f4_th8_2 = PARS_UNBN[par_index[41]]*torch.square( (PARS_UNBN[par_index[39]]-PARS_UNBN[par_index[43]]-TH8) )    #low
    f4_th8_3 = PARS_UNBN[par_index[41]]*torch.square( (PARS_UNBN[par_index[39]]+PARS_UNBN[par_index[43]]-TH8) )   #high
    f4_th8_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[39]]+PARS_UNBN[par_index[43]],f4_th8_3,f4_th8_4)
    f4_th8 = torch.where(TH8<=PARS_UNBN[par_index[39]]+PARS_UNBN[par_index[42]],f4_th8_1,f4_th8)
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[39]]-PARS_UNBN[par_index[42]],f4_th8_2,f4_th8)
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[39]]-PARS_UNBN[par_index[43]],f4_th8_4,f4_th8)
    
    energy_hydr = f1*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    
    energy_hydr = torch.sum(energy_hydr, dim=1)
    energy_hydr = energy_hydr/48
    
    
    #CROSS STACKING
    #33
    
    f2_1 = PARS_UNBN[par_index[77]]*0.5*(torch.square( HYDR_R-PARS_UNBN[par_index[78]])-torch.square( PARS_UNBN[par_index[79]]-PARS_UNBN[par_index[78]]))
    f2_2 = PARS_UNBN[par_index[77]]*PARS_UNBN[par_index[80]]*torch.square( (HYDR_R-PARS_UNBN[par_index[84]]) )  #blow
    f2_3 = PARS_UNBN[par_index[77]]*PARS_UNBN[par_index[81]]*torch.square( (HYDR_R-PARS_UNBN[par_index[85]]) )   #bhigh
    f2_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f2 = torch.where(HYDR_R<PARS_UNBN[par_index[85]],f2_3,f2_4)
    f2 = torch.where(HYDR_R<=PARS_UNBN[par_index[83]],f2_1,f2)
    f2 = torch.where(HYDR_R<PARS_UNBN[par_index[82]],f2_2,f2)
    f2 = torch.where(HYDR_R<PARS_UNBN[par_index[84]],f2_4,f2)
   
    f4_th1_1 = 1.-PARS_UNBN[par_index[87]]*torch.square( (TH1-PARS_UNBN[par_index[86]]) )
    f4_th1_2 = PARS_UNBN[par_index[88]]*torch.square( (PARS_UNBN[par_index[86]]-PARS_UNBN[par_index[90]]-TH1) )    #low
    f4_th1_3 = PARS_UNBN[par_index[88]]*torch.square( (PARS_UNBN[par_index[86]]+PARS_UNBN[par_index[90]]-TH1) )   #high
    f4_th1_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[86]]+PARS_UNBN[par_index[90]],f4_th1_3,f4_th1_4)
    f4_th1 = torch.where(TH1<=PARS_UNBN[par_index[86]]+PARS_UNBN[par_index[89]],f4_th1_1,f4_th1)
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[86]]-PARS_UNBN[par_index[89]],f4_th1_2,f4_th1)
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[86]]-PARS_UNBN[par_index[90]],f4_th1_4,f4_th1)
   
    f4_th2_1 = 1.-PARS_UNBN[par_index[92]]*torch.square( (TH2-PARS_UNBN[par_index[91]]) )
    f4_th2_2 = PARS_UNBN[par_index[93]]*torch.square( (PARS_UNBN[par_index[91]]-PARS_UNBN[par_index[95]]-TH2) )    #low
    f4_th2_3 = PARS_UNBN[par_index[93]]*torch.square( (PARS_UNBN[par_index[91]]+PARS_UNBN[par_index[95]]-TH2) )   #high
    f4_th2_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[91]]+PARS_UNBN[par_index[95]],f4_th2_3,f4_th2_4)
    f4_th2 = torch.where(TH2<=PARS_UNBN[par_index[91]]+PARS_UNBN[par_index[94]],f4_th2_1,f4_th2)
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[91]]-PARS_UNBN[par_index[94]],f4_th2_2,f4_th2)
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[91]]-PARS_UNBN[par_index[95]],f4_th2_4,f4_th2)
   
    f4_th3_1 = 1.-PARS_UNBN[par_index[97]]*torch.square( (TH3-PARS_UNBN[par_index[96]]) )
    f4_th3_2 = PARS_UNBN[par_index[98]]*torch.square( (PARS_UNBN[par_index[96]]-PARS_UNBN[par_index[100]]-TH3) )    #low
    f4_th3_3 = PARS_UNBN[par_index[98]]*torch.square( (PARS_UNBN[par_index[96]]+PARS_UNBN[par_index[100]]-TH3) )   #high
    f4_th3_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[96]]+PARS_UNBN[par_index[100]],f4_th3_3,f4_th3_4)
    f4_th3 = torch.where(TH3<=PARS_UNBN[par_index[96]]+PARS_UNBN[par_index[99]],f4_th3_1,f4_th3)
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[96]]-PARS_UNBN[par_index[99]],f4_th3_2,f4_th3)
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[96]]-PARS_UNBN[par_index[100]],f4_th3_4,f4_th3)
   
    f4_th4_1 = 1.-PARS_UNBN[par_index[102]]*torch.square( (TH4_UNBN-PARS_UNBN[par_index[101]]) )
    f4_th4_2 = PARS_UNBN[par_index[103]]*torch.square( (PARS_UNBN[par_index[101]]-PARS_UNBN[par_index[105]]-TH4_UNBN) )    #low
    f4_th4_3 = PARS_UNBN[par_index[103]]*torch.square( (PARS_UNBN[par_index[101]]+PARS_UNBN[par_index[105]]-TH4_UNBN) )   #high
    f4_th4_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[101]]+PARS_UNBN[par_index[105]],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_UNBN<=PARS_UNBN[par_index[101]]+PARS_UNBN[par_index[104]],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[101]]-PARS_UNBN[par_index[104]],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[101]]-PARS_UNBN[par_index[105]],f4_th4_4,f4_th4)
   
    f4_th7_1 = 1.-PARS_UNBN[par_index[107]]*torch.square( (TH7-PARS_UNBN[par_index[106]]) )
    f4_th7_2 = PARS_UNBN[par_index[108]]*torch.square( (PARS_UNBN[par_index[106]]-PARS_UNBN[par_index[110]]-TH7) )    #low
    f4_th7_3 = PARS_UNBN[par_index[108]]*torch.square( (PARS_UNBN[par_index[106]]+PARS_UNBN[par_index[110]]-TH7) )   #high
    f4_th7_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[106]]+PARS_UNBN[par_index[110]],f4_th7_3,f4_th7_4)
    f4_th7 = torch.where(TH7<=PARS_UNBN[par_index[106]]+PARS_UNBN[par_index[109]],f4_th7_1,f4_th7)
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[106]]-PARS_UNBN[par_index[109]],f4_th7_2,f4_th7)
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[106]]-PARS_UNBN[par_index[110]],f4_th7_4,f4_th7)
   
    f4_th8_1 = 1.-PARS_UNBN[par_index[112]]*torch.square( (TH8-PARS_UNBN[par_index[111]]) )
    f4_th8_2 = PARS_UNBN[par_index[113]]*torch.square( (PARS_UNBN[par_index[111]]-PARS_UNBN[par_index[115]]-TH8) )    #low
    f4_th8_3 = PARS_UNBN[par_index[113]]*torch.square( (PARS_UNBN[par_index[111]]+PARS_UNBN[par_index[115]]-TH8) )   #high
    f4_th8_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
   
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[111]]+PARS_UNBN[par_index[115]],f4_th8_3,f4_th8_4)
    f4_th8 = torch.where(TH8<=PARS_UNBN[par_index[111]]+PARS_UNBN[par_index[114]],f4_th8_1,f4_th8)
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[111]]-PARS_UNBN[par_index[114]],f4_th8_2,f4_th8)
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[111]]-PARS_UNBN[par_index[115]],f4_th8_4,f4_th8)
    
    energy_crst = f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    
    #55
    f2_1 = PARS_UNBN[par_index[116]]*0.5*(torch.square( HYDR_R-PARS_UNBN[par_index[117]])-torch.square( PARS_UNBN[par_index[118]]-PARS_UNBN[par_index[117]]))
    f2_2 = PARS_UNBN[par_index[116]]*PARS_UNBN[par_index[119]]*torch.square( (HYDR_R-PARS_UNBN[par_index[123]]) )  #blow
    f2_3 = PARS_UNBN[par_index[116]]*PARS_UNBN[par_index[120]]*torch.square( (HYDR_R-PARS_UNBN[par_index[124]]) )   #bhigh
    f2_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f2 = torch.where(HYDR_R<PARS_UNBN[par_index[124]],f2_3,f2_4)
    f2 = torch.where(HYDR_R<=PARS_UNBN[par_index[122]],f2_1,f2)
    f2 = torch.where(HYDR_R<PARS_UNBN[par_index[121]],f2_2,f2)
    f2 = torch.where(HYDR_R<PARS_UNBN[par_index[123]],f2_4,f2)
    
    f4_th1_1 = 1.-PARS_UNBN[par_index[126]]*torch.square( (TH1-PARS_UNBN[par_index[125]]) )
    f4_th1_2 = PARS_UNBN[par_index[127]]*torch.square( (PARS_UNBN[par_index[125]]-PARS_UNBN[par_index[129]]-TH1) )    #low
    f4_th1_3 = PARS_UNBN[par_index[127]]*torch.square( (PARS_UNBN[par_index[125]]+PARS_UNBN[par_index[129]]-TH1) )   #high
    f4_th1_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[125]]+PARS_UNBN[par_index[129]],f4_th1_3,f4_th1_4)
    f4_th1 = torch.where(TH1<=PARS_UNBN[par_index[125]]+PARS_UNBN[par_index[128]],f4_th1_1,f4_th1)
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[125]]-PARS_UNBN[par_index[128]],f4_th1_2,f4_th1)
    f4_th1 = torch.where(TH1<PARS_UNBN[par_index[125]]-PARS_UNBN[par_index[129]],f4_th1_4,f4_th1)
    
    f4_th2_1 = 1.-PARS_UNBN[par_index[131]]*torch.square( (TH2-PARS_UNBN[par_index[130]]) )
    f4_th2_2 = PARS_UNBN[par_index[132]]*torch.square( (PARS_UNBN[par_index[130]]-PARS_UNBN[par_index[134]]-TH2) )    #low
    f4_th2_3 = PARS_UNBN[par_index[132]]*torch.square( (PARS_UNBN[par_index[130]]+PARS_UNBN[par_index[134]]-TH2) )   #high
    f4_th2_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[130]]+PARS_UNBN[par_index[134]],f4_th2_3,f4_th2_4)
    f4_th2 = torch.where(TH2<=PARS_UNBN[par_index[130]]+PARS_UNBN[par_index[133]],f4_th2_1,f4_th2)
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[130]]-PARS_UNBN[par_index[133]],f4_th2_2,f4_th2)
    f4_th2 = torch.where(TH2<PARS_UNBN[par_index[130]]-PARS_UNBN[par_index[134]],f4_th2_4,f4_th2)
    
    f4_th3_1 = 1.-PARS_UNBN[par_index[136]]*torch.square( (TH3-PARS_UNBN[par_index[135]]) )
    f4_th3_2 = PARS_UNBN[par_index[137]]*torch.square( (PARS_UNBN[par_index[135]]-PARS_UNBN[par_index[139]]-TH3) )    #low
    f4_th3_3 = PARS_UNBN[par_index[137]]*torch.square( (PARS_UNBN[par_index[135]]+PARS_UNBN[par_index[139]]-TH3) )   #high
    f4_th3_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[135]]+PARS_UNBN[par_index[139]],f4_th3_3,f4_th3_4)
    f4_th3 = torch.where(TH3<=PARS_UNBN[par_index[135]]+PARS_UNBN[par_index[138]],f4_th3_1,f4_th3)
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[135]]-PARS_UNBN[par_index[138]],f4_th3_2,f4_th3)
    f4_th3 = torch.where(TH3<PARS_UNBN[par_index[135]]-PARS_UNBN[par_index[139]],f4_th3_4,f4_th3)
    
    f4_th4_1 = 1.-PARS_UNBN[par_index[141]]*torch.square( (TH4_UNBN-PARS_UNBN[par_index[140]]) )
    f4_th4_2 = PARS_UNBN[par_index[142]]*torch.square( (PARS_UNBN[par_index[140]]-PARS_UNBN[par_index[144]]-TH4_UNBN) )    #low
    f4_th4_3 = PARS_UNBN[par_index[142]]*torch.square( (PARS_UNBN[par_index[140]]+PARS_UNBN[par_index[144]]-TH4_UNBN) )   #high
    f4_th4_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[140]]+PARS_UNBN[par_index[144]],f4_th4_3,f4_th4_4)
    f4_th4 = torch.where(TH4_UNBN<=PARS_UNBN[par_index[140]]+PARS_UNBN[par_index[143]],f4_th4_1,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[140]]-PARS_UNBN[par_index[143]],f4_th4_2,f4_th4)
    f4_th4 = torch.where(TH4_UNBN<PARS_UNBN[par_index[140]]-PARS_UNBN[par_index[144]],f4_th4_4,f4_th4)
    
    f4_th7_1 = 1.-PARS_UNBN[par_index[146]]*torch.square( (TH7-PARS_UNBN[par_index[145]]) )
    f4_th7_2 = PARS_UNBN[par_index[147]]*torch.square( (PARS_UNBN[par_index[145]]-PARS_UNBN[par_index[149]]-TH7) )    #low
    f4_th7_3 = PARS_UNBN[par_index[147]]*torch.square( (PARS_UNBN[par_index[145]]+PARS_UNBN[par_index[149]]-TH7) )   #high
    f4_th7_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)
    
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[145]]+PARS_UNBN[par_index[149]],f4_th7_3,f4_th7_4)
    f4_th7 = torch.where(TH7<=PARS_UNBN[par_index[145]]+PARS_UNBN[par_index[148]],f4_th7_1,f4_th7)
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[145]]-PARS_UNBN[par_index[148]],f4_th7_2,f4_th7)
    f4_th7 = torch.where(TH7<PARS_UNBN[par_index[145]]-PARS_UNBN[par_index[149]],f4_th7_4,f4_th7)
    
    f4_th8_1 = 1.-PARS_UNBN[par_index[151]]*torch.square( (TH8-PARS_UNBN[par_index[150]]) )
    f4_th8_2 = PARS_UNBN[par_index[152]]*torch.square( (PARS_UNBN[par_index[150]]-PARS_UNBN[par_index[154]]-TH8) )    #low
    f4_th8_3 = PARS_UNBN[par_index[152]]*torch.square( (PARS_UNBN[par_index[150]]+PARS_UNBN[par_index[154]]-TH8) )   #high
    f4_th8_4 = torch.zeros(len(types_unbn),len(types_unbn[0]),device=device)

    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[150]]+PARS_UNBN[par_index[154]],f4_th8_3,f4_th8_4)
    f4_th8 = torch.where(TH8<=PARS_UNBN[par_index[150]]+PARS_UNBN[par_index[153]],f4_th8_1,f4_th8)
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[150]]-PARS_UNBN[par_index[153]],f4_th8_2,f4_th8)
    f4_th8 = torch.where(TH8<PARS_UNBN[par_index[150]]-PARS_UNBN[par_index[154]],f4_th8_4,f4_th8)
    
    energy_crst = energy_crst +  f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    
    energy_crst = torch.sum(energy_crst, dim=1)
    energy_crst = energy_crst/48

    energy_tot = energy_fene+energy_stck+energy_hydr+energy_crst

#print(counts)

print("Energy_fene in:", energy_fene.device)

test_file = open("out_test.txt", 'w')
print(energy_fene, file=test_file)

print(energy_stck, file=test_file)

print(energy_hydr, file=test_file)

print(energy_crst, file=test_file)

print(energy_tot, file=test_file)

#print(f1, file=test_file)
#print(f4_th4, file=test_file)
#print(f4_th5, file=test_file)
#print(f4_th6, file=test_file)
#print(f5_phi1, file=test_file)
#print(f5_phi2, file=test_file)

test_file.close()