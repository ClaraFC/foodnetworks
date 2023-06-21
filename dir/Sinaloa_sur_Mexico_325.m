food_web_filename = 'Sinaloa_sur_Mexico_325';
crc = 'crc_6bf2c7f7';
n = 37;
l = 36;
b0 = [0.0208772887 1.37911928 1.303709 1.39769471 0.6781125 2.60085964 1.38055933 0.346753269 1.02232158 0.5638489 2.79353857 1.02099669 0.9595986 0.7183409 0.7257744 0.5608707 0.839779854 2.00379419 10.4765673 1.028559 5.826627 0.710176349 2.959223 3.84638143 6.67976046 8.227674 4.28952169 3.72297931 1.20495093 0.8140599 1.00495 6.30602026 6.75974751 25.3970032 124.624 109.321259 4.1185].';
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3193.0164288 998.52580674489 1e-30].';
q = [0 0.084 0.0158 0.0155 0.0729 0.0083 0.2221 3e-04 0.0038 0.0063 0.0109 0.0054 0.0209 0.0135 0.001 0.0122 7e-04 0.214 4.7681 0.253635615 0.0164 0.00308 0.30651 0.004778 0.2049 0.0139 0.00044 0.009606641 0.886999965 0.000753 3.0464 0.000265 0.00468 0 0 0 1695.92632961795].';
r = [1.45363235 6.39974165 7.29731464 1.77589464 1.2318126 6.61240864 4.78486252 1.727355 6.873436 3.27565265 9.876839 4.062565 2.84709978 4.8323 1.90498817 2.054789 1.78620934 12.4013844 97.84571 12.4356775 31.5213814 4.040352 33.12789 34.664505 41.97495 13.1593037 8.333692 19.3527069 5.058386 3.7207408 11.4198895 88.12491 27.12599 1273.68958 532.1694048 166.420967790815 0].';
F = [
0 0 0 0 0 0 0 0 0.2676903 0 0 0 0 0 0 0 0.0001868527 0 0.7515114 0 0.06002304 0 0.2514035 0 0 0 0 0 0 0 0 0 0 0.4968152 0 0 0
0 0 0.3714372 0.8281761 0 0.1978529 0 0 0.3524023 0 0 0.01770962 0 0 0 0 0 0.7039682 2.999126 0.3053733 0.3813553 0 0 0 0 0 0 0 0.3120966 0.08100336 0.2156525 2.341534 0.3829335 0 0 0 0
0 0 0.1709101 0 0.2992936 0 0 0 0 0 0.1487658 0.02862596 0 0 0 0 0 0 0 0 1.064655 0 0 0.3926495 0.6632643 0 0 0 0.2638485 0.1190535 0.06639486 8.553498 0.1459317 0.2110897 0 0.4232585 0
0 0.1537896 0 0 0 0 0 0 0.1007316 0 0 0 0 0 0 0 0 0.4737017 1.584519 0.2740863 0.01941015 0 0 0 0.02704637 0 0 0 0 0.05816415 0.1954495 0.309718 0.1691147 0.05894843 0 0.03123098 0
0 0 0 0.034234 0 0 0 0 0 0.009261745 0 0 0 0.02447464 0.0434036 0 0.0002072147 0.04092236 0.4403164 0 0.04253247 0 0.02163147 0 0.007825828 0.7772014 0 0.00678203 0.1374393 0.05422777 0.05262979 0.2026437 0.2148727 0.04198517 0 0.01769284 0
0 0.602282 0.0835547 0 0 1.336847 0 0 0 0 0 0 0 0 0 0 0 0 3.109948 0 5.34573 0 1.529675 0 0.02689871 0 0 0 0 0.2574444 0.2121108 0.06571728 0.4612962 0.07034537 0 0 0
0 0 0 0 0 0 0 0.1484953 0.7355837 0 0 0.03832528 0.07016668 0 0 0.8441611 0 0 5.019908 0 0.0007262621 0.6590563 0.01162867 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0.01518659 0.0482976 0 0 0 0.07837481 0.0283823 0.2414929 0.0002502877 0.00497016 0 0.05024585 0 0.04250099 0 0 0.01536299 0.08154849 0.03718948 0.2417486 0 1.613849 0 0 0 0 0 0.001025196 0.01637209 0.001054537 0 0.002523307 0 0 0 0
0 0 0 0 0 0 0.161826 0.04853269 0 0.007666041 0 0 0.2717931 0.7829179 0.1455203 0.01493467 0.001345998 1.351561 3.313187 0.1351931 0.05551272 0 0 1.162332 0.1830995 0 0 0 0.04972752 0 0.146642 0 0.09320875 3.991059 0.007679152 0.08192716 0.2627167
0 0 0 0 0 0 0 0 0 0.054833 0 0 0 0 0 0 0 0 0 0 1.340024 0 0 0 0 0 0 0 0 0.4322071 0.2528321 0.1685095 1.714044 0.9702342 0 0 0
0 0 0 0 0 0 0.3112702 0 0 0 0.1293927 0.669831 0.2878216 0 0.7877513 0 0.003261821 3.121038 3.725989 0 0.6561129 0 2.53576 0.2675329 0.211979 0.03825783 0 0 0.420785 0.006581969 0.1731075 1.69427 0.133197 0.4760579 0 0.6483021 0
0 0 0 0 0 0.5051709 0 0 0 0 0 0.130562 0 0 0 0 0 0 0 0 0.1337395 0 0.01396498 0.02461986 0.09318429 0 0 0 0.5851911 0.3025038 0.2570727 0.9339344 0.510625 1.724201 0 1.317182 0
0.006482944 0 0 0 0 0 0 0 0 0 0 0 0.2011385 0 0 0 0.002426398 0 0 0 0.6824771 0 0 0 0.1061672 0.9018791 0 0 0 0.08956101 0.0384683 1.081208 0.227238 0.3827794 1.568721 0.6537 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0.4821688 0 0 0 0 0 0 0.6089494 0 0 1.156283 0 0 0 0 0 0.1921302 0.2395594 1.98483 0.3797192 2.883166 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.1826596 0 1.148478 0.8524786 0 0.3219267 0 0 0 0 0.5050119 0.1829295 0.5417969 0 0.6417508 0
0 0 0 0 0 1.163238 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.388977 0 0.1303238 0 0 0 0 0 0 0.01995726 0.0233259 0 0.05425047 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3816135 0 0 0 0 0 0.1432924 0.0002456024 0 0.1029675 0.01465609 0 0 0 0.0002531498 0.01463907 0.0002223175 0 0.003684659 1.602601 0 0.07753631 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3376285 0.3833177 1.855042 0.3802392 0 0 0 0 3.452077 0.2336489 0.9207496 0 0 17.04584
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4991837 1.145481 118.303 47.99091 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3181979 0.02940587 0 0 0 0 0 0 2.88104 0.3280058 0.7477256 0 1.99284 10.55885
0 0 1.221419 0 0 0 0.2408195 0 0.6990604 0 0 0 0.8467184 0 0 0 0.03839088 0 3.556297 0 3.050496 0 0 0 0.9481533 0.4918116 0 3.814367 0.06622496 0 0.1965133 0.8623692 0.2809947 17.29744 0 28.96851 0.5563839
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3299388 0 0 0 0 0 0 0 0 0 0 0 6.555114 0 0 0.9592017
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.222821 0 0 0.6463451 2.36981 0 0 0 0 0 0 0 0.1633058 0 0 46.1308 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 34.18651 26.2462
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20.7163 39.32988 0 12.81544
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4369709 0.08946382 1.55787 5.091627 0.201042 0 0 0.03921007 0.4243597 0.0386312 0.287673 0.07633006 14.4952 3.457356
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4867515 0.2127962 1.67124 17.94195
0 0 0 0 0 0 0 0 0 0 1.635179 0 0 0 0 0 0 0 0 0 0 0 0 4.715353 6.069303 1.390662 0 0 0 0 0 0 0.3778628 0 0 15.32959 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3900757 0 0 0.5076657 0.6127717 0 0 0 0.4688693 0 0.1150872 1.5123 0.06625348 0 0 2.239081 5.092254
0 0 0 0 0 0 0 0 0 0 0 0.0923609 0 0 0 0 0 0 0 0 0.100492 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6.800911
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4356653 0.01708263 0 2.93782 2.449807 15.32256
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.9409279 0.04607345 0 0 17.9085 139.3405
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.2109426 0 8.947757 0 35.53349
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.2765254 181.0422 1691.072 0 324.6112
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.367515 2.235621 3.342505 0.8021058 0.5662781 3.278033 1.728799 0.5673102 2.994184 1.57872 4.492255 1.486571 1.346363 1.787796 1.070254 0.853988 0.5089828 6.286606 36.51962 3.414905 14.89186 1.982577 11.62064 15.19275 20.36276 6.009809 6.184813 6.133443 2.753513 1.628424 4.308011 41.26202 9.86445 517.3735 868.7032 708.971 0
];

