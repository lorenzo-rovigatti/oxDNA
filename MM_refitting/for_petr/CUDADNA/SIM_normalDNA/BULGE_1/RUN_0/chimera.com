set bg_color white
~bond #0
bond #0:0
bond #0:1
bond #0:2
bond #0:3
bond #0:4
bond #0:5
bond #0:6
bond #0:7
bond #0:8
bond #0:9
bond #0:10
bond #0:11
bond #0:12
bond #0:13
bond #0:14
bond #0:15
bond #0:16
bond #0:17
bond #0:18
bond #0:19
bond #0:20
bond #0:0,1@C
bond #0:1,2@C
bond #0:2,3@C
bond #0:3,4@C
bond #0:4,5@C
bond #0:5,6@C
bond #0:6,7@C
bond #0:7,8@C
bond #0:8,9@C
bond #0:9,10@C
bond #0:11,12@C
bond #0:12,13@C
bond #0:13,14@C
bond #0:14,15@C
bond #0:15,16@C
bond #0:16,17@C
bond #0:17,18@C
bond #0:18,19@C
bond #0:19,20@C
color sandy brown #0:ALA
col deep sky blue #0:ALA@P
col deep sky blue #0:ALA@O
col deep sky blue #0:ALA@S
col deep sky blue #0:ALA@K
bondcolor sandy brown #0:ALA
color blue #0:GLY
col deep sky blue #0:GLY@P
col deep sky blue #0:GLY@O
col deep sky blue #0:GLY@S
col deep sky blue #0:GLY@K
bondcolor blue #0:GLY
color red #0:CYS
col deep sky blue #0:CYS@P
col deep sky blue #0:CYS@O
col deep sky blue #0:CYS@S
col deep sky blue #0:CYS@K
bondcolor red #0:CYS
color green #0:TYR
col deep sky blue #0:TYR@P
col deep sky blue #0:TYR@O
col deep sky blue #0:TYR@S
col deep sky blue #0:TYR@K
bondcolor green #0:TYR
color purple #0:ARG
col deep sky blue #0:ARG@P
col deep sky blue #0:ARG@O
col deep sky blue #0:ARG@S
col deep sky blue #0:ARG@K
bondcolor purple #0:ARG
color dark gray #0:PHE
col deep sky blue #0:PHE@P
col deep sky blue #0:PHE@O
col deep sky blue #0:PHE@S
col deep sky blue #0:PHE@K
bondcolor dark gray #0:PHE
color yellow #0:LYS
col deep sky blue #0:LYS@P
col deep sky blue #0:LYS@O
col deep sky blue #0:LYS@S
col deep sky blue #0:LYS@K
bondcolor yellow #0:LYS
color orange #0:SER
col deep sky blue #0:SER@P
col deep sky blue #0:SER@O
col deep sky blue #0:SER@S
col deep sky blue #0:SER@K
bondcolor orange #0:SER
color deep pink #0:PRO
col deep sky blue #0:PRO@P
col deep sky blue #0:PRO@O
col deep sky blue #0:PRO@S
col deep sky blue #0:PRO@K
bondcolor deep pink #0:PRO
color magenta #0:VAL
col deep sky blue #0:VAL@P
col deep sky blue #0:VAL@O
col deep sky blue #0:VAL@S
col deep sky blue #0:VAL@K
bondcolor magenta #0:VAL
color sienna #0:ASN
col deep sky blue #0:ASN@P
col deep sky blue #0:ASN@O
col deep sky blue #0:ASN@S
col deep sky blue #0:ASN@K
bondcolor sienna #0:ASN
color goldenrod #0:ASP
col deep sky blue #0:ASP@P
col deep sky blue #0:ASP@O
col deep sky blue #0:ASP@S
col deep sky blue #0:ASP@K
bondcolor goldenrod #0:ASP
color gray #0:CYX
col deep sky blue #0:CYX@P
col deep sky blue #0:CYX@O
col deep sky blue #0:CYX@S
col deep sky blue #0:CYX@K
bondcolor gray #0:CYX
color plum #0:HSP
col deep sky blue #0:HSP@P
col deep sky blue #0:HSP@O
col deep sky blue #0:HSP@S
col deep sky blue #0:HSP@K
bondcolor plum #0:HSP
color olive drab #0:HSD
col deep sky blue #0:HSD@P
col deep sky blue #0:HSD@O
col deep sky blue #0:HSD@S
col deep sky blue #0:HSD@K
bondcolor olive drab #0:HSD
color dark red #0:MET
col deep sky blue #0:MET@P
col deep sky blue #0:MET@O
col deep sky blue #0:MET@S
col deep sky blue #0:MET@K
bondcolor dark red #0:MET
color steel blue #0:LEU
col deep sky blue #0:LEU@P
col deep sky blue #0:LEU@O
col deep sky blue #0:LEU@S
col deep sky blue #0:LEU@K
bondcolor steel blue #0:LEU
aniso scale 0.75 smoothing 4
setattr m stickScale 0.6 #0
