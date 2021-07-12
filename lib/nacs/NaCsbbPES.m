function V = NaCsbbPES(r)

wavenum2hartree = 4.55634011e-6;
bohr2angstrom = 0.529177210903;
h = 6.62607003912593e-34;
hartree = 4.35974464404181e-18;
% r = r*bohr2angstrom;

% cs_d2 = 11732.3071049*wavenum2hartree;
% cs_d1 = 11732.3071049*wavenum2hartree;
cs_d1 = h*(335116048808294)/hartree;

R = [4.5
    4.6
    4.7
    4.8
    4.9
    5
    5.1
    5.2
    5.3
    5.4
    5.5
    5.6
    5.7
    5.8
    5.9
    6
    6.1
    6.2
    6.3
    6.4
    6.5
    6.6
    6.7
    6.8
    6.9
    7
    7.1
    7.2
    7.3
    7.4
    7.5
    7.6
    7.7
    7.8
    7.9
    8
    8.1
    8.2
    8.3
    8.4
    8.5
    8.6
    8.7
    8.8
    8.9
    9
    9.1
    9.2
    9.3
    9.4
    9.5
    9.6
    9.7
    9.8
    9.9
    10
    10.1
    10.2
    10.3
    10.4
    10.5
    10.6
    10.7
    10.8
    10.9
    11
    11.1
    11.2
    11.3
    11.4
    11.5
    11.6
    11.7
    11.8
    11.9
    12
    12.1
    12.2
    12.3
    12.4
    12.5
    12.6
    12.7
    12.8
    12.9
    13
    13.1
    13.2
    13.3
    13.4
    13.5
    13.6
    13.7
    13.8
    13.9
    14
    14.1
    14.2
    14.3
    14.4
    14.5
    14.6
    14.7
    14.8
    14.9
    15
    15.1
    15.2
    15.3
    15.4
    15.5
    15.6
    15.7
    15.8
    15.9
    16
    16.1
    16.2
    16.3
    16.4
    16.5
    16.6
    16.7
    16.8
    16.9
    17
    17.1
    17.2
    17.3
    17.4
    17.5
    17.6
    17.7
    17.8
    17.9
    18
    18.1
    18.2
    18.3
    18.4
    18.5
    18.6
    18.7
    18.8
    18.9
    19
    19.1
    19.2
    19.3
    19.4
    19.5
    19.6
    19.7
    19.8
    19.9
    20
    20.1
    20.2
    20.3
    20.4
    20.5
    20.6
    20.7
    20.8
    20.9
    21
    21.1
    21.2
    21.3
    21.4
    21.5
    21.6
    21.7
    21.8
    21.9
    22
    22.1
    22.2
    22.3
    22.4
    22.5
    22.6
    22.7
    22.8
    22.9
    23
    23.1
    23.2
    23.3
    23.4
    23.5
    23.6
    23.7
    23.8
    23.9
    24
    24.1
    24.2
    24.3
    24.4
    24.5
    24.6
    24.7
    24.8
    24.9
    25
    25.1
    25.2
    25.3
    25.4
    25.5
    25.6
    25.7
    25.8
    25.9
    26
    26.1
    26.2
    26.3
    26.4
    26.5
    26.6
    26.7
    26.8
    26.9
    27
    27.1
    27.2
    27.3
    27.4
    27.5
    27.6
    27.7
    27.8
    27.9
    28
    28.1
    28.2
    28.3
    28.4
    28.5
    28.6
    28.7
    28.8
    28.9
    29
    29.1
    29.2
    29.3
    29.4
    29.5
    29.6
    29.7
    29.8
    29.9
    30
    30.1
    30.2
    30.3
    30.4
    30.5
    30.6
    30.7
    30.8
    30.9
    31
    31.1
    31.2
    31.3
    31.4
    31.5
    31.6
    31.7
    31.8
    31.9
    32
    32.1
    32.2
    32.3
    32.4
    32.5
    32.6
    32.7
    32.8
    32.9
    33
    33.1
    33.2
    33.3
    33.4
    33.5
    33.6
    33.7
    33.8
    33.9
    34
    34.1
    34.2
    34.3
    34.4
    34.5
    34.6
    34.7
    34.8
    34.9
    35
    35.1
    35.2
    35.3
    35.4
    35.5
    35.6
    35.7
    35.8
    35.9
    36
    36.1
    36.2
    36.3
    36.4
    36.5
    36.6
    36.7
    36.8
    36.9
    37
    37.1
    37.2
    37.3
    37.4
    37.5
    37.6
    37.7
    37.8
    37.9
    38
    38.1
    38.2
    38.3
    38.4
    38.5
    38.6
    38.7
    38.8
    38.9
    39
    39.1
    39.2
    39.3
    39.4
    39.5
    39.6
    39.7
    39.8
    39.9
    40
    40.1
    40.2
    40.3
    40.4
    40.5
    40.6
    40.7
    40.8
    40.9
    41
    41.1
    41.2
    41.3
    41.4
    41.5
    41.6
    41.7
    41.8
    41.9
    42
    42.1
    42.2
    42.3
    42.4
    42.5
    42.6
    42.7
    42.8
    42.9
    43
    43.1
    43.2
    43.3
    43.4
    43.5
    43.6
    43.7
    43.8
    43.9
    44
    44.1
    44.2
    44.3
    44.4
    44.5
    44.6
    44.7
    44.8
    44.9
    45
    45.1
    45.2
    45.3
    45.4
    45.5
    45.6
    45.7
    45.8
    45.9
    46
    46.1
    46.2
    46.3
    46.4
    46.5
    46.6
    46.7
    46.8
    46.9
    47
    47.1
    47.2
    47.3
    47.4
    47.5
    47.6
    47.7
    47.8
    47.9
    48
    48.1
    48.2
    48.3
    48.4
    48.5
    48.6
    48.7
    48.8
    48.9
    49
    49.1
    49.2
    49.3
    49.4
    49.5
    49.6
    49.7
    49.8
    49.9
    50
    50.1
    50.2
    50.3
    50.4
    50.5
    50.6
    50.7
    50.8
    50.9
    51
    51.1
    51.2
    51.3
    51.4
    51.5
    51.6
    51.7
    51.8
    51.9
    52
    52.1
    52.2
    52.3
    52.4
    52.5
    52.6
    52.7
    52.8
    52.9
    53
    53.1
    53.2
    53.3
    53.4
    53.5
    53.6
    53.7
    53.8
    53.9
    54
    54.1
    54.2
    54.3
    54.4
    54.5
    54.6
    54.7
    54.8
    54.9
    55
    55.1
    55.2
    55.3
    55.4
    55.5
    55.6
    55.7
    55.8
    55.9
    56
    56.1
    56.2
    56.3
    56.4
    56.5
    56.6
    56.7
    56.8
    56.9
    57
    57.1
    57.2
    57.3
    57.4
    57.5
    57.6
    57.7
    57.8
    57.9
    58
    58.1
    58.2
    58.3
    58.4
    58.5
    58.6
    58.7
    58.8
    58.9
    59
    59.1
    59.2
    59.3
    59.4
    59.5
    59.6
    59.7
    59.8
    59.9
    60
    60.1
    60.2
    60.3
    60.4
    60.5
    60.6
    60.7
    60.8
    60.9
    61
    61.1
    61.2
    61.3
    61.4
    61.5
    61.6
    61.7
    61.8
    61.9
    62
    62.1
    62.2
    62.3
    62.4
    62.5
    62.6
    62.7
    62.8
    62.9
    63
    63.1
    63.2
    63.3
    63.4
    63.5
    63.6
    63.7
    63.8
    63.9
    64
    64.1
    64.2
    64.3
    64.4
    64.5
    64.6
    64.7
    64.8
    64.9
    65
    65.1
    65.2
    65.3
    65.4
    65.5
    65.6
    65.7
    65.8
    65.9
    66
    66.1
    66.2
    66.3
    66.4
    66.5
    66.6
    66.7
    66.8
    66.9
    67
    67.1
    67.2
    67.3
    67.4
    67.5
    67.6
    67.7
    67.8
    67.9
    68
    68.1
    68.2
    68.3
    68.4
    68.5
    68.6
    68.7
    68.8
    68.9
    69
    69.1
    69.2
    69.3
    69.4
    69.5
    69.6
    69.7
    69.8
    69.9
    70
    70.1
    70.2
    70.3
    70.4
    70.5
    70.6
    70.7
    70.8
    70.9
    71
    71.1
    71.2
    71.3
    71.4
    71.5
    71.6
    71.7
    71.8
    71.9
    72
    72.1
    72.2
    72.3
    72.4
    72.5
    72.6
    72.7
    72.8
    72.9
    73
    73.1
    73.2
    73.3
    73.4
    73.5
    73.6
    73.7
    73.8
    73.9
    74
    74.1
    74.2
    74.3
    74.4
    74.5
    74.6
    74.7
    74.8
    74.9
    75
    75.1
    75.2
    75.3
    75.4
    75.5
    75.6
    75.7
    75.8
    75.9
    76
    76.1
    76.2
    76.3
    76.4
    76.5
    76.6
    76.7
    76.8
    76.9
    77
    77.1
    77.2
    77.3
    77.4
    77.5
    77.6
    77.7
    77.8
    77.9
    78
    78.1
    78.2
    78.3
    78.4
    78.5
    78.6
    78.7
    78.8
    78.9
    79
    79.1
    79.2
    79.3
    79.4
    79.5
    79.6
    79.7
    79.8
    79.9
    80
    80.1
    80.2
    80.3
    80.4
    80.5
    80.6
    80.7
    80.8
    80.9
    81
    81.1
    81.2
    81.3
    81.4
    81.5
    81.6
    81.7
    81.8
    81.9
    82
    82.1
    82.2
    82.3
    82.4
    82.5
    82.6
    82.7
    82.8
    82.9
    83
    83.1
    83.2
    83.3
    83.4
    83.5
    83.6
    83.7
    83.8
    83.9
    84
    84.1
    84.2
    84.3
    84.4
    84.5
    84.6
    84.7
    84.8
    84.9
    85
    85.1
    85.2
    85.3
    85.4
    85.5
    85.6
    85.7
    85.8
    85.9
    86
    86.1
    86.2
    86.3
    86.4
    86.5
    86.6
    86.7
    86.8
    86.9
    87
    87.1
    87.2
    87.3
    87.4
    87.5
    87.6
    87.7
    87.8
    87.9
    88
    88.1
    88.2
    88.3
    88.4
    88.5
    88.6
    88.7
    88.8
    88.9
    89
    89.1
    89.2
    89.3
    89.4
    89.5
    89.6
    89.7
    89.8
    89.9
    90
    90.1
    90.2
    90.3
    90.4
    90.5
    90.6
    90.7
    90.8
    90.9
    91
    91.1
    91.2
    91.3
    91.4
    91.5
    91.6
    91.7
    91.8
    91.9
    92
    92.1
    92.2
    92.3
    92.4
    92.5
    92.6
    92.7
    92.8
    92.9
    93
    93.1
    93.2
    93.3
    93.4
    93.5
    93.6
    93.7
    93.8
    93.9
    94
    94.1
    94.2
    94.3
    94.4
    94.5
    94.6
    94.7
    94.8
    94.9
    95
    95.1
    95.2
    95.3
    95.4
    95.5
    95.6
    95.7
    95.8
    95.9
    96
    96.1
    96.2
    96.3
    96.4
    96.5
    96.6
    96.7
    96.8
    96.9
    97
    97.1
    97.2
    97.3
    97.4
    97.5
    97.6
    97.7
    97.8
    97.9
    98
    98.1
    98.2
    98.3
    98.4
    98.5
    98.6
    98.7
    98.8
    98.9
    99
    99.1
    99.2
    99.3
    99.4
    99.5
    99.6
    99.7
    99.8
    99.9
    100
    105
    110
    115
    120
    125
    130
    135
    140
    145
    150
    155
    160
    165
    170
    175
    180
    185
    190
    195
    200
    205
    210
    215
    220
    225
    230
    235
    240
    245
    250
    255
    260
    265
    270
    275
    280
    285
    290
    295
    300
    305
    310
    315
    320
    325
    330
    335
    340
    345
    350
    355
    360
    365
    370
    375
    380
    385
    390
    395
    400
    405
    410
    415
    420
    425
    430
    435
    440
    445
    450
    455
    460
    465
    470
    475
    480
    485
    490
    495
    500
    550
    600
    650
    700
    750
    800
    850
    900
    950
    1000
    1050
    1100
    1150
    1200
    1250
    1300
    1350
    1400
    1450
    1500
    1550
    1600
    1650
    1700
    1750
    1800
    1850
    1900
    1950
    2000
    2050
    2100
    2150
    2200
    2250
    2300
    2350
    2400
    2450
    2500
    2550
    2600
    2650
    2700
    2750
    2800
    2850
    2900
    2950
    3000
    3050
    3100
    3150
    3200
    3250
    3300
    3350
    3400
    3450
    3500
    3550
    3600
    3650
    3700
    3750
    3800
    3850
    3900
    3950
    4000
    4050
    4100
    4150
    4200
    4250
    4300
    4350
    4400
    4450
    4500
    4550
    4600
    4650
    4700
    4750
    4800
    4850
    4900
    4950
    5000
    5500
    6000
    6500
    7000
    7500
    8000
    8500
    9000
    9500
    10000
    20000
    30000
    40000
    50000
    60000
    70000
    80000
    90000
    100000];

U = [0.037745789
    0.029818028
    0.02282662
    0.016686805
    0.011313828
    0.006622929
    0.002529353
    -0.001051659
    -0.00420486
    -0.007006126
    -0.009503983
    -0.01174171
    -0.013761377
    -0.015588526
    -0.017236735
    -0.018719124
    -0.020043796
    -0.021213489
    -0.022230744
    -0.023102413
    -0.023842087
    -0.024463928
    -0.02497907
    -0.025391493
    -0.025704224
    -0.025924433
    -0.026065495
    -0.026141
    -0.02615796
    -0.026116896
    -0.026017965
    -0.025860948
    -0.025645686
    -0.025376376
    -0.025061968
    -0.024711409
    -0.024330373
    -0.023922217
    -0.023490251
    -0.023037304
    -0.022564353
    -0.022071918
    -0.021560524
    -0.021030721
    -0.020483958
    -0.019922719
    -0.019349548
    -0.01876699
    -0.018177619
    -0.017584205
    -0.016989584
    -0.016396594
    -0.015808069
    -0.015226201
    -0.014651537
    -0.014084362
    -0.013524961
    -0.01297363
    -0.012430965
    -0.01189793
    -0.011375509
    -0.010864687
    -0.010366421
    -0.009881385
    -0.009410089
    -0.008953037
    -0.008510738
    -0.008083655
    -0.007672068
    -0.007276212
    -0.006896319
    -0.00653262
    -0.006185131
    -0.00585346
    -0.005537166
    -0.005235809
    -0.004948952
    -0.004676255
    -0.004417463
    -0.004172323
    -0.003940585
    -0.00372196
    -0.00351588
    -0.003321659
    -0.003138608
    -0.00296604
    -0.002803264
    -0.002649594
    -0.002504367
    -0.002367161
    -0.002237696
    -0.002115687
    -0.002000854
    -0.001892914
    -0.001791586
    -0.001696576
    -0.001607527
    -0.001524068
    -0.001445826
    -0.001372428
    -0.001303502
    -0.001238677
    -0.001177614
    -0.001120008
    -0.001065555
    -0.001013964
    -0.000965055
    -0.000918713
    -0.000874822
    -0.000833268
    -0.000793937
    -0.000756728
    -0.000721561
    -0.000688359
    -0.000657043
    -0.000627534
    -0.000599749
    -0.000573596
    -0.000548982
    -0.000525814
    -0.000503999
    -0.000483424
    -0.000463938
    -0.000445389
    -0.00042763
    -0.000410582
    -0.000394218
    -0.000378509
    -0.000363425
    -0.000348938
    -0.000335019
    -0.000321644
    -0.000308797
    -0.000296468
    -0.000284643
    -0.000273312
    -0.000262463
    -0.000252083
    -0.000242161
    -0.000232686
    -0.000223645
    -0.000215027
    -0.00020682
    -0.00019901
    -0.000191581
    -0.000184515
    -0.000177797
    -0.000171409
    -0.000165335
    -0.000159558
    -0.000154062
    -0.000148829
    -0.000143843
    -0.000139088
    -0.000134546
    -0.000130202
    -0.000126037
    -0.000122037
    -0.000118188
    -0.000114487
    -0.00011093
    -0.00010751
    -0.000104224
    -0.000101066
    -9.80336E-05
    -9.51205E-05
    -9.23225E-05
    -8.96349E-05
    -8.70533E-05
    -8.45729E-05
    -8.21891E-05
    -7.98973E-05
    -7.76929E-05
    -7.55713E-05
    -7.35278E-05
    -7.15578E-05
    -6.96568E-05
    -6.782E-05
    -6.60433E-05
    -6.43245E-05
    -6.26618E-05
    -6.10535E-05
    -5.94981E-05
    -5.79937E-05
    -5.65387E-05
    -5.51314E-05
    -5.37701E-05
    -5.24531E-05
    -5.11787E-05
    -4.99453E-05
    -4.87511E-05
    -4.75945E-05
    -4.64737E-05
    -4.53871E-05
    -4.43329E-05
    -4.33096E-05
    -4.23153E-05
    -4.13485E-05
    -4.04073E-05
    -3.94902E-05
    -3.85954E-05
    -3.77213E-05
    -3.68666E-05
    -3.6031E-05
    -3.52141E-05
    -3.44157E-05
    -3.36354E-05
    -3.28729E-05
    -3.21278E-05
    -3.13998E-05
    -3.06887E-05
    -2.9994E-05
    -2.93154E-05
    -2.86527E-05
    -2.80055E-05
    -2.73735E-05
    -2.67563E-05
    -2.61536E-05
    -2.55652E-05
    -2.49906E-05
    -2.44295E-05
    -2.38817E-05
    -2.33468E-05
    -2.28244E-05
    -2.23143E-05
    -2.18161E-05
    -2.13295E-05
    -2.08542E-05
    -2.03899E-05
    -1.99366E-05
    -1.9494E-05
    -1.90619E-05
    -1.86401E-05
    -1.82283E-05
    -1.78265E-05
    -1.74345E-05
    -1.70519E-05
    -1.66787E-05
    -1.63146E-05
    -1.59594E-05
    -1.56129E-05
    -1.5275E-05
    -1.49455E-05
    -1.46241E-05
    -1.43106E-05
    -1.40049E-05
    -1.37068E-05
    -1.3416E-05
    -1.31324E-05
    -1.28557E-05
    -1.25859E-05
    -1.23227E-05
    -1.2066E-05
    -1.18157E-05
    -1.15715E-05
    -1.13333E-05
    -1.1101E-05
    -1.08743E-05
    -1.06532E-05
    -1.04375E-05
    -1.02269E-05
    -1.00214E-05
    -9.82E-06
    -9.62E-06
    -9.43E-06
    -9.25E-06
    -9.06E-06
    -8.89E-06
    -8.71E-06
    -8.54E-06
    -8.37E-06
    -8.21E-06
    -8.05E-06
    -7.89E-06
    -7.74E-06
    -7.59E-06
    -7.45E-06
    -7.30E-06
    -7.17E-06
    -7.03E-06
    -6.90E-06
    -6.77E-06
    -6.64E-06
    -6.52E-06
    -6.39E-06
    -6.28E-06
    -6.16E-06
    -6.05E-06
    -5.94E-06
    -5.83E-06
    -5.72E-06
    -5.62E-06
    -5.52E-06
    -5.42E-06
    -5.32E-06
    -5.23E-06
    -5.14E-06
    -5.04E-06
    -4.96E-06
    -4.87E-06
    -4.78E-06
    -4.70E-06
    -4.62E-06
    -4.54E-06
    -4.46E-06
    -4.38E-06
    -4.30E-06
    -4.23E-06
    -4.16E-06
    -4.08E-06
    -4.01E-06
    -3.94E-06
    -3.88E-06
    -3.81E-06
    -3.74E-06
    -3.68E-06
    -3.61E-06
    -3.55E-06
    -3.49E-06
    -3.43E-06
    -3.37E-06
    -3.32E-06
    -3.26E-06
    -3.20E-06
    -3.15E-06
    -3.10E-06
    -3.04E-06
    -2.99E-06
    -2.94E-06
    -2.89E-06
    -2.84E-06
    -2.80E-06
    -2.75E-06
    -2.70E-06
    -2.66E-06
    -2.61E-06
    -2.57E-06
    -2.53E-06
    -2.49E-06
    -2.45E-06
    -2.41E-06
    -2.37E-06
    -2.33E-06
    -2.29E-06
    -2.26E-06
    -2.22E-06
    -2.19E-06
    -2.15E-06
    -2.12E-06
    -2.09E-06
    -2.06E-06
    -2.03E-06
    -2.00E-06
    -1.97E-06
    -1.94E-06
    -1.91E-06
    -1.88E-06
    -1.85E-06
    -1.83E-06
    -1.80E-06
    -1.78E-06
    -1.75E-06
    -1.73E-06
    -1.71E-06
    -1.68E-06
    -1.66E-06
    -1.64E-06
    -1.62E-06
    -1.60E-06
    -1.57E-06
    -1.55E-06
    -1.53E-06
    -1.51E-06
    -1.49E-06
    -1.48E-06
    -1.46E-06
    -1.44E-06
    -1.42E-06
    -1.40E-06
    -1.38E-06
    -1.37E-06
    -1.35E-06
    -1.33E-06
    -1.31E-06
    -1.30E-06
    -1.28E-06
    -1.26E-06
    -1.25E-06
    -1.23E-06
    -1.22E-06
    -1.20E-06
    -1.18E-06
    -1.17E-06
    -1.15E-06
    -1.14E-06
    -1.12E-06
    -1.11E-06
    -1.09E-06
    -1.08E-06
    -1.06E-06
    -1.05E-06
    -1.04E-06
    -1.02E-06
    -1.01E-06
    -9.95E-07
    -9.82E-07
    -9.69E-07
    -9.56E-07
    -9.43E-07
    -9.31E-07
    -9.19E-07
    -9.07E-07
    -8.95E-07
    -8.83E-07
    -8.72E-07
    -8.61E-07
    -8.50E-07
    -8.39E-07
    -8.28E-07
    -8.17E-07
    -8.07E-07
    -7.97E-07
    -7.87E-07
    -7.77E-07
    -7.67E-07
    -7.58E-07
    -7.48E-07
    -7.39E-07
    -7.30E-07
    -7.21E-07
    -7.12E-07
    -7.04E-07
    -6.95E-07
    -6.87E-07
    -6.78E-07
    -6.70E-07
    -6.62E-07
    -6.54E-07
    -6.46E-07
    -6.39E-07
    -6.31E-07
    -6.24E-07
    -6.16E-07
    -6.09E-07
    -6.02E-07
    -5.95E-07
    -5.88E-07
    -5.81E-07
    -5.74E-07
    -5.67E-07
    -5.60E-07
    -5.54E-07
    -5.47E-07
    -5.41E-07
    -5.34E-07
    -5.28E-07
    -5.22E-07
    -5.16E-07
    -5.10E-07
    -5.04E-07
    -4.98E-07
    -4.92E-07
    -4.86E-07
    -4.81E-07
    -4.75E-07
    -4.70E-07
    -4.64E-07
    -4.59E-07
    -4.53E-07
    -4.48E-07
    -4.43E-07
    -4.38E-07
    -4.33E-07
    -4.28E-07
    -4.23E-07
    -4.18E-07
    -4.13E-07
    -4.09E-07
    -4.04E-07
    -3.99E-07
    -3.95E-07
    -3.90E-07
    -3.86E-07
    -3.82E-07
    -3.77E-07
    -3.73E-07
    -3.69E-07
    -3.65E-07
    -3.61E-07
    -3.57E-07
    -3.53E-07
    -3.49E-07
    -3.45E-07
    -3.41E-07
    -3.37E-07
    -3.34E-07
    -3.30E-07
    -3.26E-07
    -3.23E-07
    -3.19E-07
    -3.16E-07
    -3.12E-07
    -3.09E-07
    -3.06E-07
    -3.02E-07
    -2.99E-07
    -2.96E-07
    -2.93E-07
    -2.89E-07
    -2.86E-07
    -2.83E-07
    -2.80E-07
    -2.77E-07
    -2.74E-07
    -2.71E-07
    -2.69E-07
    -2.66E-07
    -2.63E-07
    -2.60E-07
    -2.57E-07
    -2.55E-07
    -2.52E-07
    -2.49E-07
    -2.47E-07
    -2.44E-07
    -2.42E-07
    -2.39E-07
    -2.37E-07
    -2.34E-07
    -2.32E-07
    -2.30E-07
    -2.27E-07
    -2.25E-07
    -2.23E-07
    -2.20E-07
    -2.18E-07
    -2.16E-07
    -2.14E-07
    -2.12E-07
    -2.09E-07
    -2.07E-07
    -2.05E-07
    -2.03E-07
    -2.01E-07
    -1.99E-07
    -1.97E-07
    -1.95E-07
    -1.93E-07
    -1.91E-07
    -1.89E-07
    -1.87E-07
    -1.85E-07
    -1.84E-07
    -1.82E-07
    -1.80E-07
    -1.78E-07
    -1.76E-07
    -1.75E-07
    -1.73E-07
    -1.71E-07
    -1.70E-07
    -1.68E-07
    -1.66E-07
    -1.65E-07
    -1.63E-07
    -1.61E-07
    -1.60E-07
    -1.58E-07
    -1.57E-07
    -1.55E-07
    -1.54E-07
    -1.52E-07
    -1.51E-07
    -1.49E-07
    -1.48E-07
    -1.46E-07
    -1.45E-07
    -1.44E-07
    -1.42E-07
    -1.41E-07
    -1.40E-07
    -1.38E-07
    -1.37E-07
    -1.36E-07
    -1.34E-07
    -1.33E-07
    -1.32E-07
    -1.30E-07
    -1.29E-07
    -1.28E-07
    -1.27E-07
    -1.26E-07
    -1.24E-07
    -1.23E-07
    -1.22E-07
    -1.21E-07
    -1.20E-07
    -1.19E-07
    -1.18E-07
    -1.17E-07
    -1.15E-07
    -1.14E-07
    -1.13E-07
    -1.12E-07
    -1.11E-07
    -1.10E-07
    -1.09E-07
    -1.08E-07
    -1.07E-07
    -1.06E-07
    -1.05E-07
    -1.04E-07
    -1.03E-07
    -1.02E-07
    -1.01E-07
    -1.01E-07
    -9.96E-08
    -9.87E-08
    -9.79E-08
    -9.70E-08
    -9.61E-08
    -9.52E-08
    -9.44E-08
    -9.35E-08
    -9.27E-08
    -9.19E-08
    -9.11E-08
    -9.03E-08
    -8.95E-08
    -8.87E-08
    -8.79E-08
    -8.71E-08
    -8.63E-08
    -8.56E-08
    -8.48E-08
    -8.41E-08
    -8.34E-08
    -8.26E-08
    -8.19E-08
    -8.12E-08
    -8.05E-08
    -7.98E-08
    -7.91E-08
    -7.84E-08
    -7.77E-08
    -7.71E-08
    -7.64E-08
    -7.57E-08
    -7.51E-08
    -7.45E-08
    -7.38E-08
    -7.32E-08
    -7.26E-08
    -7.19E-08
    -7.13E-08
    -7.07E-08
    -7.01E-08
    -6.95E-08
    -6.89E-08
    -6.83E-08
    -6.78E-08
    -6.72E-08
    -6.66E-08
    -6.60E-08
    -6.55E-08
    -6.49E-08
    -6.44E-08
    -6.38E-08
    -6.33E-08
    -6.28E-08
    -6.23E-08
    -6.17E-08
    -6.12E-08
    -6.07E-08
    -6.02E-08
    -5.97E-08
    -5.92E-08
    -5.87E-08
    -5.82E-08
    -5.77E-08
    -5.73E-08
    -5.68E-08
    -5.63E-08
    -5.59E-08
    -5.54E-08
    -5.49E-08
    -5.45E-08
    -5.40E-08
    -5.36E-08
    -5.32E-08
    -5.27E-08
    -5.23E-08
    -5.19E-08
    -5.15E-08
    -5.10E-08
    -5.06E-08
    -5.02E-08
    -4.98E-08
    -4.94E-08
    -4.90E-08
    -4.86E-08
    -4.82E-08
    -4.78E-08
    -4.75E-08
    -4.71E-08
    -4.67E-08
    -4.63E-08
    -4.60E-08
    -4.56E-08
    -4.52E-08
    -4.49E-08
    -4.45E-08
    -4.42E-08
    -4.38E-08
    -4.35E-08
    -4.31E-08
    -4.28E-08
    -4.25E-08
    -4.21E-08
    -4.18E-08
    -4.15E-08
    -4.11E-08
    -4.08E-08
    -4.05E-08
    -4.02E-08
    -3.99E-08
    -3.96E-08
    -3.93E-08
    -3.90E-08
    -3.86E-08
    -3.83E-08
    -3.81E-08
    -3.78E-08
    -3.75E-08
    -3.72E-08
    -3.69E-08
    -3.66E-08
    -3.63E-08
    -3.60E-08
    -3.58E-08
    -3.55E-08
    -3.52E-08
    -3.50E-08
    -3.47E-08
    -3.44E-08
    -3.42E-08
    -3.39E-08
    -3.36E-08
    -3.34E-08
    -3.31E-08
    -3.29E-08
    -3.26E-08
    -3.24E-08
    -3.21E-08
    -3.19E-08
    -3.17E-08
    -3.14E-08
    -3.12E-08
    -3.10E-08
    -3.07E-08
    -3.05E-08
    -3.03E-08
    -3.01E-08
    -2.98E-08
    -2.96E-08
    -2.94E-08
    -2.92E-08
    -2.90E-08
    -2.87E-08
    -2.85E-08
    -2.83E-08
    -2.81E-08
    -2.79E-08
    -2.77E-08
    -2.75E-08
    -2.73E-08
    -2.71E-08
    -2.69E-08
    -2.67E-08
    -2.65E-08
    -2.63E-08
    -2.61E-08
    -2.59E-08
    -2.58E-08
    -2.56E-08
    -2.54E-08
    -2.52E-08
    -2.50E-08
    -2.48E-08
    -2.47E-08
    -2.45E-08
    -2.43E-08
    -2.41E-08
    -2.40E-08
    -2.38E-08
    -2.36E-08
    -2.35E-08
    -2.33E-08
    -2.31E-08
    -2.30E-08
    -2.28E-08
    -2.26E-08
    -2.25E-08
    -2.23E-08
    -2.22E-08
    -2.20E-08
    -2.18E-08
    -2.17E-08
    -2.15E-08
    -2.14E-08
    -2.12E-08
    -2.11E-08
    -2.09E-08
    -2.08E-08
    -2.06E-08
    -2.05E-08
    -2.04E-08
    -2.02E-08
    -2.01E-08
    -1.99E-08
    -1.98E-08
    -1.97E-08
    -1.95E-08
    -1.94E-08
    -1.93E-08
    -1.91E-08
    -1.90E-08
    -1.89E-08
    -1.87E-08
    -1.86E-08
    -1.85E-08
    -1.83E-08
    -1.82E-08
    -1.81E-08
    -1.80E-08
    -1.78E-08
    -1.77E-08
    -1.76E-08
    -1.75E-08
    -1.74E-08
    -1.73E-08
    -1.71E-08
    -1.70E-08
    -1.69E-08
    -1.68E-08
    -1.67E-08
    -1.66E-08
    -1.65E-08
    -1.63E-08
    -1.62E-08
    -1.61E-08
    -1.60E-08
    -1.59E-08
    -1.58E-08
    -1.57E-08
    -1.56E-08
    -1.55E-08
    -1.54E-08
    -1.53E-08
    -1.52E-08
    -1.51E-08
    -1.50E-08
    -1.49E-08
    -1.48E-08
    -1.47E-08
    -1.46E-08
    -1.45E-08
    -1.44E-08
    -1.43E-08
    -1.42E-08
    -1.41E-08
    -1.40E-08
    -1.39E-08
    -1.38E-08
    -1.38E-08
    -1.37E-08
    -1.36E-08
    -1.35E-08
    -1.34E-08
    -1.33E-08
    -1.32E-08
    -1.31E-08
    -1.31E-08
    -1.30E-08
    -1.29E-08
    -1.28E-08
    -1.27E-08
    -1.26E-08
    -1.26E-08
    -1.25E-08
    -1.24E-08
    -1.23E-08
    -1.22E-08
    -1.22E-08
    -1.21E-08
    -1.20E-08
    -1.19E-08
    -1.18E-08
    -1.18E-08
    -1.17E-08
    -1.16E-08
    -1.15E-08
    -1.15E-08
    -1.14E-08
    -1.13E-08
    -1.13E-08
    -1.12E-08
    -1.11E-08
    -1.10E-08
    -1.10E-08
    -1.09E-08
    -1.08E-08
    -1.08E-08
    -1.07E-08
    -1.06E-08
    -1.06E-08
    -1.05E-08
    -1.04E-08
    -1.04E-08
    -1.03E-08
    -1.02E-08
    -1.02E-08
    -1.01E-08
    -1.01E-08
    -9.99E-09
    -9.93E-09
    -9.87E-09
    -9.81E-09
    -9.75E-09
    -9.69E-09
    -9.63E-09
    -9.57E-09
    -9.51E-09
    -9.45E-09
    -9.39E-09
    -9.34E-09
    -9.28E-09
    -9.22E-09
    -9.17E-09
    -9.11E-09
    -9.05E-09
    -9.00E-09
    -8.95E-09
    -8.89E-09
    -8.84E-09
    -8.78E-09
    -8.73E-09
    -8.68E-09
    -8.63E-09
    -8.57E-09
    -8.52E-09
    -8.47E-09
    -8.42E-09
    -8.37E-09
    -8.32E-09
    -8.27E-09
    -6.16E-09
    -4.66E-09
    -3.56E-09
    -2.76E-09
    -2.15E-09
    -1.70E-09
    -1.36E-09
    -1.09E-09
    -8.81E-10
    -7.18E-10
    -5.89E-10
    -4.87E-10
    -4.04E-10
    -3.38E-10
    -2.84E-10
    -2.40E-10
    -2.03E-10
    -1.73E-10
    -1.48E-10
    -1.27E-10
    -1.09E-10
    -9.47E-11
    -8.22E-11
    -7.16E-11
    -6.25E-11
    -5.48E-11
    -4.81E-11
    -4.24E-11
    -3.74E-11
    -3.32E-11
    -2.94E-11
    -2.62E-11
    -2.33E-11
    -2.09E-11
    -1.87E-11
    -1.68E-11
    -1.51E-11
    -1.36E-11
    -1.22E-11
    -1.11E-11
    -1.00E-11
    -9.09E-12
    -8.25E-12
    -7.51E-12
    -6.84E-12
    -6.24E-12
    -5.70E-12
    -5.22E-12
    -4.78E-12
    -4.38E-12
    -4.02E-12
    -3.70E-12
    -3.40E-12
    -3.14E-12
    -2.89E-12
    -2.67E-12
    -2.47E-12
    -2.29E-12
    -2.12E-12
    -1.96E-12
    -1.82E-12
    -1.69E-12
    -1.57E-12
    -1.46E-12
    -1.36E-12
    -1.27E-12
    -1.19E-12
    -1.11E-12
    -1.03E-12
    -9.67E-13
    -9.05E-13
    -8.47E-13
    -7.94E-13
    -7.45E-13
    -6.99E-13
    -6.56E-13
    -6.17E-13
    -5.80E-13
    -5.45E-13
    -5.13E-13
    -2.90E-13
    -1.72E-13
    -1.06E-13
    -6.80E-14
    -4.49E-14
    -3.05E-14
    -2.12E-14
    -1.50E-14
    -1.09E-14
    -7.99E-15
    -5.96E-15
    -4.51E-15
    -3.45E-15
    -2.67E-15
    -2.09E-15
    -1.65E-15
    -1.32E-15
    -1.06E-15
    -8.58E-16
    -7.00E-16
    -5.75E-16
    -4.75E-16
    -3.95E-16
    -3.30E-16
    -2.78E-16
    -2.34E-16
    -1.99E-16
    -1.69E-16
    -1.45E-16
    -1.25E-16
    -1.07E-16
    -9.29E-17
    -8.07E-17
    -7.03E-17
    -6.14E-17
    -5.38E-17
    -4.73E-17
    -4.17E-17
    -3.69E-17
    -3.26E-17
    -2.90E-17
    -2.58E-17
    -2.30E-17
    -2.06E-17
    -1.84E-17
    -1.65E-17
    -1.49E-17
    -1.34E-17
    -1.21E-17
    -1.09E-17
    -9.90E-18
    -8.98E-18
    -8.16E-18
    -7.42E-18
    -6.76E-18
    -6.17E-18
    -5.64E-18
    -5.16E-18
    -4.72E-18
    -4.33E-18
    -3.98E-18
    -3.66E-18
    -3.37E-18
    -3.10E-18
    -2.86E-18
    -2.65E-18
    -2.45E-18
    -2.26E-18
    -2.10E-18
    -1.94E-18
    -1.80E-18
    -1.68E-18
    -1.56E-18
    -1.45E-18
    -1.35E-18
    -1.26E-18
    -1.18E-18
    -1.10E-18
    -1.03E-18
    -9.59E-19
    -8.98E-19
    -8.41E-19
    -7.88E-19
    -7.39E-19
    -6.93E-19
    -6.51E-19
    -6.12E-19
    -5.75E-19
    -5.41E-19
    -5.10E-19
    -2.88E-19
    -1.71E-19
    -1.06E-19
    -6.77E-20
    -4.47E-20
    -3.04E-20
    -2.11E-20
    -1.50E-20
    -1.08E-20
    -7.96E-21
    -1.24E-22
    -1.09E-23
    -1.94E-24
    -5.09E-25
    -1.71E-25
    -6.77E-26
    -3.04E-26
    -1.50E-26
    -7.96E-27];


U1 =@(r) interp1(R,U,r,'spline');

V = U1(r) + cs_d1;

end