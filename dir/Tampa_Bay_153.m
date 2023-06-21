food_web_filename = 'Tampa_Bay_153';
crc = 'crc_0ff48146';
n = 52;
l = 51;
b0 = [0.000216968474 0.0184970144 0.227183774 0.0984339043 0.0200000033 0.000273548416 0.00415257 0.0272239055 0.108115189 0.300000072 9.09441e-05 0.0259648785 0.220000058 1.96918627e-05 0.00252293865 0.100000076 0.06484298 0.522366643 2.79999924 1.24145038e-06 0.0183000043 0.009787013 0.08899991 0.023921689 0.3850111 0.32 0.69 0.09 0.08 0.1815 0.0389 0.120790623 0.13644 0.0565 0.037 1.25845969 0.08 0.0824824944 0.122 0.9 0.05 0.1 0.0256353635 4 6.48041 1.20710933 0.7861352 20 29.778 175.617 25 151.303154349018].';
p = [0 0 0 0 0 0 0 0 0 0.134344565982883 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 893.34 1899.6139656 5463.9 0].';
q = [0 0 0.02 0.035 0 0 0 0 0.005 5e-04 0 0 0.05 0 0 0.0165 0 0.25 0.5 0 0.005 0 0.01 0 0.02 0.001 0 0.001 0.001 0 0 0.01 0 0.01 0 0 0.001 0 0.099 0 0 0 0 0 0 0 0 0 0 0 0 6604.80192239974].';
r = [0.00334344758 0.0557570755 0.224959552 0.0569512472 0.008799996 0.0016456017 0.00510196434 0.0274836235 0.0661310554 0.07019998 0.00114017026 0.04698044 0.127599925 0.000498127833 0.0146948807 0.145999938 2.50509644 6.675989 15.6799955 7.702851e-05 0.07868999 0.112283833 0.284800172 0.023921689 3.338046 1.72191989 5.865 0.5220001 0.646783948 1.83315 0.338040978 1.20790625 1.43261993 0.29832 0.3108 10.5710611 0.5056 0.210330352 0.381250024 4.95 0.275 0.24 0.115359135 2 56.3795624 21.969389 13.6787529 120 148.89 316.6023276 910.65 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.002214145 0.002214145 0.0005535362 0 0 0 0 0.0005535362
0.0001159389 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0005796944 0 0 0 0 0 0 0 0 0.0005796944 0 0 0 0.001159389 0 0 0.001159389 0 0 0.001159389 0 0 0 0.01739083 0.005796944 0 0 0 0.04625961 0.01078232 0.007999783 0 0 0 0 0.0229559
0 0 0 0 0 5.367221e-06 0 0 0 0 0 0 0 0 0 0 0.0005367275 0 0 0 0 0.0005367275 0 0 0.0005367275 0.02689005 0.01073455 0.002683638 0.0005367275 0.05533661 0.0005367275 0 0.02689005 0 0 0.1095998 0.05533661 0 0 0.08050913 0.02683638 0.02683638 0 0 0.08501764 0.02742678 0 0 0 0 0 0
0 0.0002949802 0 0 0 1.474901e-05 0 0 0 0 0 0.003657754 0 0 0 0 0.0002949802 0 0 0 0 0 0 0 0.01461627 0.03362774 0.004527946 0.004527946 0.01461627 0.004527946 0.0002949802 0.003215284 0 0 0 0.004527946 0.02942427 0 0 0.01461627 0.007315508 0.007374504 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 4.796115e-05 0 0 0.0002170873 9.59223e-05 0.001201553 0.0002170873 0.00098699 0 0 0.0007572813 9.59223e-05 0 0.00500563 0.0005048542 9.59223e-05 0.001153592 9.59223e-05 9.59223e-05 0.00500563 0.01001378 9.59223e-05 2.524271e-05 9.59223e-05 0 9.59223e-05 0 9.59223e-05 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000575674 0.002307494 0.001731819 0 5.756741e-05 0 4.797284e-06 0.0001151348
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0003190843 0 0 0 0 0 0.001963596 0.001325427 0 0 0 0 0 0 0 0.008958906 0.01195339 0 0 0 0 0 2.454495e-05
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0004307245 0.0009332363 0 0 0 0 0 0 0.001363961 0.0004307245 0.002297197 0.009253398 0 0.004659003 0 0.004659003 0.004659003 0 0 0.0139124 0.01979179 0 0.004659003 0.00358937 0.0007178741 0 0 0 0.0004307245 0 0 0 0 0 0 0
0 0 0 0 0 0.0001637666 0 0 0 0 0 0 0 0 0 0 0.0009825996 0.001801433 0 0 0 0 0 0 0.00931832 0.004585465 0.004585465 0.0009825996 0 0.004585465 0 0.02458137 0.00833572 0 0 0.02782394 0.04046673 0 0.001801433 0.01310133 0.00818833 0.00818833 0.003275332 0 0.0009825996 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0009701031 0.01346385 0 0 0 0 0 0 0.01184702 0.006820118 0.00470353 0.001616838 0 0.0004997501 0 0.03371843 0.0004997501 0.0001763824 0 0.02869153 0.0480054 0 0.001793221 0 0 0 0.005879412 0 0.0009701031 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.105188e-06 0 0.0004210377 0 0 0 0 0 0 0 0 0.0004210377 0.0008420754 0.0004210377 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.005187372 0 0 0.005187372 0.005187372 0 0 0 0.003114506 0 0 0.0425302 0.03319709 0 0 0.005187372 0.001041641 0 0.002083282 0 0 0.001447881 0 0 0 0 0 0
0 0 0 0 0 3.520351e-05 0 0 0 0 0 0.001056106 0 0 0 0 0.04224422 0.0003520352 0 0 0 0 0 0 0.05498789 0.03559075 0.002816281 0.01056105 0 0.0003520352 0 0.0003520352 0.04481407 0 0 0.03520352 0.0701254 0.0003520352 0.0003520352 0.01760176 0 0.01760176 0.01760176 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.491319e-06 0 0 0 0 0 0 0 0 0.0001313852 0.0003500871 0.0002627703 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 2.215522e-05 0 0 0 0 0 0 0 0 0 0 0.009446988 0 0.0002215522 0.0001107761 0.0001107761 0.0001107761 0 0 0 0 0 0.0005538806 0.0005538806 0 0 0 0 0 0 0 0.01102222 0 0 0 0 0 0 0
0 0 0 0 0 2.69973e-05 0 0 0 0 0 0 0 0 0.001079892 0 0 0 0 0 0 0 0 0 0.1161964 0.0005399459 0.01349865 0.03479952 0.02346065 0.006479351 0 0 0 0 0 0.03317968 0.03317968 0 0.0005399459 0 0 0 0 0 0.007019297 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.837215 1.469851 0 0 0 0 0.3673648
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.141811 2.065679 2.065679 0 1.038066 0 0 0.2090768
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.469979 0.04442588 0.04442588 0 12.32 3.114488 0 6.406679
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9.891014e-08 0 0 0 0 0 0 0 1.710047e-05 9.891014e-08 2.088103e-07 2.088103e-07 2.088103e-07 2.088103e-07 9.891014e-08 0 2.088103e-07 0 0 1.668285e-05 2.088103e-07 0 0 0 0 0 0 0 1.668285e-05 4.965289e-05 8.242512e-07 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0.0001561402 0.0008699238 0 0 0 0 1.115287e-05 1.115287e-05 0 0 0 0 0.008420416 0.04803541 0.0002899746 0.02440248 0.01262505 0.004929568 0.0001561402 0.001316039 0 0 0 0.001751001 0 0.006535582 0 0.0002899746 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01746093 0 0 0 0 0.0001794546 0 0 0 0 0 0.000699873 0 0 0 0.02618243 0.00331991 0 0 0 0.02480063 0.0523828 0.04958331 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0005446201 0 0 0 0 0 0 0 0.1329418 0 0.00108924 0.0005446201 0.001579398 0.007461295 0 0 0 0 0 0.01443243 0.003757878 0 0 0 0 0.01601183 0 0 0.2376177 0.1180192 0 0 0 0 0 0
0 0 0 0 0 0 4.784817e-05 0 0 0 0 0 0 0 0 0 0 0 0 4.784768e-07 0 0.00181823 0.00181823 0.0001913927 0.01796699 0.009100722 0.003641245 0.00181823 0.00181823 0.00181823 0.001435445 0 0.0007177225 0 0.003641245 0.0001913927 0.0001913927 0 0 0 0 0.001435445 0 0 0.0001913927 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.02695078 0 0 0 0 0 0 0 0 0.4419927 2.856782 0.4419927 0 0 0 0 1.622437
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.7560435 0.002560256 0 0 0 0 0 0 0 0.6047324 0.0460846 0.007680767 0 0.3791739 0 0.007680767 0.7560435
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3312 0 0 0 0 0 0 0 0 3.16296 1.17576 0.6624 0 0.97704 0 0.02484 1.9458
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0430067 0.01457854 0.001619838 0.001619838 0.01457854 0.01457854 0 0 0.001619838 0 0 0.1435177 0.001619838 0 0 0 0 0 0 0 0.1435177 0.2862254 0.1435177 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3521745 0.4402182 0.1760873 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.178 0.408375 0.136125 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.41623 0.06953375 0 0 0 0.00048625 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.6764275 0.5073206 0 0 0 0.3382137 0.1691069
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.2510496 1.130815 0.5261126 0 0 0 0.02401344 0.2510496
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.293571e-05 4.293527e-06 0 0 0 0 0 0.07651143 0.0004293571 0.0008587141 0.0008587141 0.006440356 0.0004293571 0.0008587141 0.0004293571 0 0.0004293571 0.005581642 0.05929421 0.02623372 0.05929421 0.1199194 0 0 0 0 0 0.0717885 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.05439442 0.2653386 0.08844621 0 0.02255378 0 0 0.01326693
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.808882 0.02265227 1.155266 0 4.666369 0 0 13.99911
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.494592 0.09984 0.09984 0 0.016896 0 0 0.824832
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.001154755 0.001732132 0 0 0 0 0 0.2552008 0 0 0 0 0 0 0.3192897
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.001037 0.005185 0.001037 0 0 0 0 0 0.60146 0.027999 0 0 0 0 0 0.400282
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.5 0 0 0 4.5 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3 0 0 0 0.2 0 0 0
0 0 0 0 0 0 0 0 0 0 0.00020002 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.04908491 0.04908491 0.02884288 0 0.04908491 0 0 0.06136614 0 0 0.1418142 0 0 0 0 0 0 0 0 0.02052205 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03845305 0 0 0 0 0 0 0 0 0.1538122 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.895978 3.906001 1.895978 0 18.87425 0.1425548 38.61809 77.23617
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.359746 0 0 0 44.0776 29.43898
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06996603 0 0.6996603 0 69.19641 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 200 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.002075976 0.05988682 0.2918216 0.05552409 0.0172 0.002900801 0.01939514 0.04430378 0.09261916 0.2233 0.0007671033 0.05226569 0.1735301 0.0002254509 0.006161149 0.1074041 1.121549 2.577501 6.218996 2.49859e-05 0.02611 0.05921326 0.2372857 0.01530988 1.467662 0.6540271 2.290822 0.1698137 0.247342 0.7375165 0.1387614 0.3961932 0.5971477 0.1204469 0.1221302 10.26903 0.6559619 0.3010611 0.4245299 3.871822 0.1717834 0.08255176 0.04806631 2 62.85998 42.97309 39.30676 76 700.7559 1579.755 4200.963 0
];
