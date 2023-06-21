food_web_filename = 'South_East_Alaska_438';
crc = 'crc_b5d4bb64';
n = 40;
l = 39;
b0 = [7e-04 0.0174518935 0.0361 2.39263341e-06 9.079065e-05 0.00144268747 0.013 0.041 3e-05 0.005954 0.01 0.251 1.11 0.0344353877 0.9239813 1.3 0.345 0.0054057166 0.2 2.38 0.6824242 2.62234569 0.89 0.28 0.26 0.6 0.198297635 1.75313652 1.02378523 0.143090144 1.37164319 2.559978 50 2.90795016 0.226004839 40 17.153183 31.8044968 23.228 13.95].';
p = [0.000154 0 0 0.000545183953406719 0.0078493701071115 0 0 0 0 0 0 0.0602400019075997 10.8831780570068 0 0.0291931891735 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3816.539616 122.64384 0].';
q = [0 0.000460451 0.001466248 0 0 0 0 0 0 0 1.18261e-06 5.6069e-06 0.5418769 0 3.05058e-10 0 0.1856 0 0 0 0 0 0.0120504741 0.000153255 0.08696734 0.000261933 0 0.000403117 0 1.47396e-07 0.000103117 0.0154440189 0 0.0311562456 2.60004e-07 0 0 0 0 1083.86866198911].';
r = [0.00614599977 0.159964055 0.311893165 0.0005396838 0.0062259296 0.04651206 0.263769984 0.93767 0.00203369977 0.321635067 0.009099999 0.210839972 6.992999 0.0227273554 2.46703 4.485 1.31444991 0.0189884249 0.18 0.4998 0.219058156 0.173074812 0.3471 0.310800016 0.280799985 0.36 0.1588364 1.38322473 2.40999055 0.1931717 3.01761484 21.5038147 580 16.38145 1.23172641 368.8 1193.86157 636.089936 20.44064 0].';
F = [
0 0 0.000462 0 3.08e-05 0.0001925 0.001001 0.0057827 0 7.7e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0.01288473 0.007847244 0.0008830658 0.003512194 0 0.0006020903 0 0 0.006582854 0.01314564 0.01204181 0.001204181 0.001204181 0.001204181 0 0.003451985 0 0.0006020903 0.03447971 0 0 0 0 0.1010508 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.0003547997 0 0.008633459 0.009421902 0.01576887 0.0001182666 0.0003153775 0.0001182666 0.0001971109 0.0001971109 0 0.0003547997 0 0 0 0 0.005795061 0.0001576887 5.124884e-06 0.0003547997 0.02432349 0.0003547997 2.444176e-05 0.2854955 0.04222116 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 5.848272e-05 0 0.0001754482 0.01312352 5.848272e-05 0.0003041101 0.008865981 0.009357235 5.848272e-05 0.001754482 0.001169654 0.001754482 0.0003508964 8.772409e-05 0.002924136 0.0002339309 0.004678618 0.001187199 0.0008772408 0.00774896 0.0006433099 0.00206444 0 0.0001462068 0 0.000859696 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0.00049725 0 0.0016575 0.06699615 0.0003315001 0.016575 0.04972501 0.04972501 0.0016575 0.01326 0.00663 0.016575 0.009945 0.0006630001 0.01369095 0.00043095 0.02111655 0.005270851 0.007326151 0.0281775 0.0001326 0.0027846 0 0.00208845 0 0.0162435 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.1182842 0.001182842 0.1774264 0.2602253 0.02957106 0.001182842 0.01774264 0.02365685 0.006742202 0.01182842 0.008634749 0.047314 7.688475e-06 0.01182842 0.00473137 0.1596837 0.1361452 0.003193674 0.06694931 0.001064558 0 0 0.09545538 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.002295 0.000255 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0.000202436 0 0 0 0 0.1102467 0.1102871 0.0001619488 0.000202436 0.001255103 0.0001214616 0.0003643848 0.001255103 0 0.001255103 0 0 0 0.0003643848 0.0009716929 0 0.0003238976 0 0 0 0.04728906 0.1145788 0.01570904 0 0 0.0002834104
0 9.229999e-05 0.0002951 0 1.3e-05 1.3e-05 0 0.0010998 0 4.679999e-05 0 0.0020163 0.0003952 0.0003431999 0.0002639 0 0.0001105 1.3e-05 0.0002717 2.86e-05 0.0001209 0.0001209 9.099999e-05 9.879999e-05 0.0002015 0.0020787 0.0006408999 0.0002821 0.0001976 0.0003431999 0.0002106 2.6e-05 7.799999e-05 5.329999e-05 0.0029146 0 0 0 0 0.0005394999
0 0 0 0 0 0 0 0 0 0 9.036001e-05 0.004518 0.01749972 0.0024096 0.0189756 0.00304212 0.01454796 6.024e-05 0.003012 9.036001e-05 0.00207828 0.00153612 0.01237932 0.002259 0.00313248 0.00117468 0.00138552 0.00475896 0.01277088 0.00725892 0.00740952 0.0213852 0.04921608 0.01123476 0.02427672 0.01198776 0 0 0.0010542 0.00141564
0 0 0 0 0 0 0 0 0 0 0 0 0.0006183624 0 0.02100099 0.01516738 0.02333443 0 0 0 0 0 0 0 0 0 0 0.002333443 0.02100099 0 0 0 0 0.003500165 0.02333443 0.6265294 0.03500164 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.0006325148 0.0001136254 0.003030011 0.001799069 0.00353375 0.0001136254 0.0003673888 0 0 0 0 0 0 0 0 0.0002878511 0 0 0.0007347777 0.001162767 0.005082844 0 0.003749639 0.013226 0.004045065 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03741662 0.05838638 0 0 0 0 0 0 0 0 0 0 0 0 0.007812262 0 0 0.3815673 1.08097 0.08182316 0.06702098 1.85973 0.507797 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.470925 0 0 1.56975 5.434325 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.02212657 0 0 1.790281 0.3783425 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0007069392 0.0001560361 0 2.229087e-05 0 0 0 0 0 0 0 0 0 0.000550903 0.0008629753 0.0001560361 0.000550903 0.00327039 0.004219344 0 0 0.01814796 0.003200333 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.02390522 0.0029994 0.000179964 0.00029994 0.003419316 0 0.000179964 0 0 0.00014997 0 0.000179964 0 0.002039592 0.004769046 0.00014997 0.000539892 0.09925015 0.03803239 0 0 0.1095681 0.01433713 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01666 0 0 0 0 0.00833 0.05831 0.09163 0 0.01666 0.5831 0.05831 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.006936842 0.03526836 0.005001828 0.0001825485 0.00116831 0 0 0 0.003687479 0 0 0 0.0008762326 0.0001095291 0.009894127 0.0008762326 0.005220886 0.07991972 0.01303396 0.008397229 0.01752465 0.1341731 0.04070831 0 0 0.002117562
0 0 0 0 0 0 0 0 0 0 0 0 0.005769161 0 0.01127871 0.03109578 0.008653741 0.000288458 0 0 0 0.007413371 0 0 0 0 0.002105744 0.007384525 0.03057655 0 0.0002048052 0.02068244 0.03441304 0.03686494 0.001557673 0.05757622 0.0110191 0 0.01367291 0.00790375
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0217516 0.0070577 0.005785 0.0005785 0 0 0 0 0 0 0 0.0108758 0 0.0217516 0.0468585 0 0.0329745 0.0293878 0.08347755 0.003471 0.01336335 0.1282535 0.00017355 0 0 0.1727401
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0891996 0.11914 0.01036 0 0 0.01554 0.0085988 0.0004144 0 0.002072 0 0.004662 0.0171458 0.08806 0.084434 0 0 0.020202 0.017094 0.0410774 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.0083772 0 0.0234468 0.0351468 0.00936 0 0.001404 0.01872 0.00936 0.003744 0.001404 0.011232 0 0.0117 0.005616 0.045864 0.0578916 0 0 0.0057096 0.0234 0.1723644 0.0136188 0.0010296 0 0 0 0.0086112
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.07577242 0.078 0.018 0.000119988 6e-04 0.06 0.00659934 0.00119988 0 0 0 0 0.003779622 0.0149985 0.1632 0.006 0.08952 0.03653635 0 0 0.01823818 0.02741726 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.01813564 0.0002647538 0.01427023 0.01180802 0.007942614 0 0.005295076 0.005295076 0 0.002647538 0.00137672 0.006618845 0 0.001588523 0.0006089338 0.02025367 0.02687251 0.0006354092 0.004050733 0.0106431 0.03417972 0.04005725 0.001773851 0.03725086 0.003415324 0.007015976 0.002726964 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03734707 0.009452037 0.01152687 0 0.002305375 0 0 0 0 0.00461075 0 0.009682574 0.0009221499 0.03780815 0.159993 0.004149674 0.06754748 0.2250046 1.468985 0.1463913 0.03112256 0.05117932 0.01729031 0 0 0.02005676
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.007448584 0.02059314 0.01402086 0.0001971684 0 0 0.001752608 0.000438152 0 0 0 0 0 0.01708793 0.03373771 0 0.005257825 0.533231 2.261741 0.3649807 8.324888e-05 0.4123011 0.6326916 0.03768108 0.03461401 0.003943368
0 0 0 0 0 0 0 0 0 0 0 0 0.00560254 0.001609925 0.000128794 0 0 0 0 0.0009337567 0.001899712 0.001674322 0.0009337567 0 0 0.0002253895 0.00257588 0.00579573 0.008468207 0.002994461 0.04069891 0.02785171 0.09089638 0.03731807 0.0237947 0.06368864 0.004861974 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.003840217 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.003840217 0 0.3867647 0.4394305 0.004937422 2.465968 2.181792 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.6655943 10.64951 0 0 5.324754 0 0 0 16.63986
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 85 85 170 85 425
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.520224 10.08089 1.764157 0 0 0.756067 0 5.040447 5.040447
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06780823 0.0180822 0.01130137 0.0006780824 0 0 0 0 0 0 0 0 0 0 0.1918973 0 0.3465001 0.107137 0.5345549 0.4084316 0.02260274 0.5510549 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 40.56 213.616 278.512 0 143.312
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1921.156 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.001554 0.04017997 0.08009548 5.500171e-06 0.001579641 0.01176515 0.066729 0.2382975 0.0005163 0.08235495 0.003808457 0.06910241 3.856487 0.007954575 0.8634605 2.10428 0.4418549 0.007082312 0.06782924 0.1943131 0.07667035 0.06057619 0.1888876 0.1132622 0.09502193 0.1432869 0.05559274 0.4841286 0.9311327 0.06761009 1.783136 6.911941 242.5299 5.229463 0.4808253 165.403 418.5377 810.7366 12.11069 0
];

