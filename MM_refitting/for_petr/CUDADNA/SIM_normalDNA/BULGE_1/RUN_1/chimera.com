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
col cyan #0:ALA@O
bondcolor sandy brown #0:ALA
color blue #0:GLY
col cyan #0:GLY@O
bondcolor blue #0:GLY
color red #0:CYS
col cyan #0:CYS@O
bondcolor red #0:CYS
color green #0:TYR
col cyan #0:TYR@O
bondcolor green #0:TYR
color yellow #0:ARG
col cyan #0:ARG@O
bondcolor yellow #0:ARG
color plum #0:PHE
col cyan #0:PHE@O
bondcolor plum #0:PHE
color sandy brown #0:LYS
col cyan #0:LYS@O
bondcolor sandy brown #0:LYS
color sandy brown #0:SER
col cyan #0:SER@O
bondcolor sandy brown #0:SER
color blue #0:PRO
col cyan #0:PRO@O
bondcolor blue #0:PRO
color red #0:VAL
col cyan #0:VAL@O
bondcolor red #0:VAL
color green #0:ASN
col cyan #0:ASN@O
bondcolor green #0:ASN
color yellow #0:ASP
col cyan #0:ASP@O
bondcolor yellow #0:ASP
aniso scale 0.75 smoothing 4
setattr m stickScale 0.6 #0
