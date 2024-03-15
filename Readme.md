# Activities

Helper functions to compute activities for perpleX pseudosections. See discussions on the list, summarized at <https://docs.google.com/document/d/1ceNjeZAwL3o1VTOb6W0hVwPLmDfhfutHEkmhXcMkKMU/edit?usp=drive_link>

Most functions are translated from MATLAB scripts by Jesse Walters.

In the perpleX world, the functions will work only if the HSC convention is used; i.e. in the HPxxxver.dat (or any other thermodynamic file used) you must comment OFF the tag HSC_conversion, BEFORE running vertex. It should look like that:

`| HSC_conversion`

The activity of TiO2 is defined by aTiO2 = exp( -(G_ru - mu_TiO2)/ RT) with \_ru = G0 + RTln(a_ru) i.e. it is a function of the Gibbs energy for rutile, the pure TiO2-bearing phase.

Therefore, we calculate it from the thermodynamics of rutile, as defined in the HP database.

This is relatively crude - if the database changes, the values will change and the code should be modified...

For now (i.e. THERMOCALC DS6.2), the thermodynamics used is defined as follows:

In perpleX format:

``` R
ru EoS = 8 | H= -944160.0
TiO2(1)
GH = -959216.6 S0 = 50.5 V0 = 1.882
c1 = 90.4 c2 = .29E-2 c5 = -623.8
b1 = .224E-4 b5 = 457.0037 b6 = 2220000 b7 = -.19E-5 b8 = 4.24
end
```

in THERMOCALC format:

``` R
ru
1 2 1.0000 10 2.0000 0
-944.36 0.05050 1.8820
0.0904 0.000002900 0.0 -0.6238
0.0000224 2220.00 4.24 -0.00190 0
```

likewise, quartz is

``` R
q         EoS = 8 | H=  -910710.0 
SiO2(1)
GH = -923062.4  S0 = 41.43  V0 = 2.269
c1 = 92.9  c2 = -.642E-3  c3 = -714900  c5 = -716.1
b5 = 525.2346  b6 = 730000  b7 = -.82E-5  b8 = 6
transition = 1  type = 4  t1 = 847  t2 = 4.95  t3 = .1188
end
```

or

``` R
q   1  1   1.0000 10   2.0000  0
-910.71   0.04143   2.2690
0.0929   -0.000000642    -714.9   -0.7161
0.0000000   730.00     6.00   -0.00820    1    847   0.00495    0.1188
```

This is, for now, hard-coded in file activity_constants.R. Brave users can try to change the values...
