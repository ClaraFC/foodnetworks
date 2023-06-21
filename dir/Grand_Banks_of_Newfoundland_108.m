food_web_filename = 'Grand_Banks_of_Newfoundland_108';
crc = 'crc_3b586ca1';
n = 50;
l = 49;
b0 = [1e-06 0.251 1e-06 0.405428 0.034 0.000227 0.013453 0.002921 0.08 0.193 0.33959 0.306902021 0.365723 0.522461 0.32995 0.376 0.547679245 0.207668 0.017 1.193 0.00202657841 0.235 0.428415179 0.295 1.58537221 0.876 0.105 0.0936 0.007718664 6.041 2.242 1.488 0.845 0.0110414578 0.9132073 1.001 0.423423648 1.60394633 0.232 1.942 0.0031659815 1.859 112.3 10.5 42.1 7.8 20.3173828 30.3666668 35.6 181].';
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2259.888 0].';
q = [0 5.79125e-05 0 0.01787302 0.001915609 1.4e-05 0.00081 0.000176 0.001995 4e-05 0.00218386552 0.00028 0.0260272622 0.00055 0.0006433468 0.00275479443 0.0009037032 0.0180450287 0.000228 0.0135664428 0.000283835 0.001632516 0 0.008827605 3.367e-06 3.7302998e-05 2.62626e-05 4.309761e-05 0.0006161116 0.024005387 0 0 0.006986532 0.0009555554 0.00128109741 0.0002444444 0.000233303 0 0.06526465 0.000348771 0.0009993271 0.0450627469 0 0 0.0328902379 0.00241818186 0 0 0 431.90275633404].';
r = [1.34169613e-05 2.34313488 9.32e-06 4.49392557 0.2912345 0.00988585 0.585878134 0.127209544 0.118447989 0.587685049 0.466621816 0.573293 0.8054794 1.55515742 0.774161637 0.472138166 0.5736705 0.4179111 0.0229670014 1.5151099 0.00263455184 0.20977667 1.12930238 0.332140446 4.179041 2.659536 0.09555 0.08180639 0.0205339417 21.26432 11.1732321 2.539123 1.9054749 0.0249978583 1.03521168 2.41187334 1.10090148 3.47521687 0.15312 1.7478 0.009991838 11.6812124 531.583252 165.6648 189.298431 58.5 167.9638 255.079956 376.648 0].';
F = [
0 0 1.686138e-08 1.686138e-08 1.686138e-08 0 0 0 0 0 0 1.686137e-07 0 0 6.74455e-08 6.74455e-08 6.74455e-08 0 0 0 0 0 3.372275e-08 0 1.686137e-07 2.69782e-07 0 3.372275e-08 0 7.419005e-07 0 0 0 0 0 0 0 0 0 2.023365e-06 1.686137e-09 2.023365e-06 8.430688e-07 1.686137e-06 5.058413e-06 3.372275e-06 0 1.686137e-07 0 0
0 0 0 0 0 0 0 0 0 0.02957337 0 0 0 0 0 0 0 0 0.002957337 0 0 0.04436005 0.04436005 0.04436005 0.04436005 0 0.01774402 0.005914674 0.002957337 1.446138 0.1537815 0 0.1596962 0 0.1626535 0.0887201 0 0.1596962 0 0 0 0 0 0 0 0 0.307563 0.2454589 0 0
0 0 0 0 0 0 0 0 1.178821e-06 1.072727e-06 0 8.251748e-08 0 1.178821e-08 8.251748e-08 3.536464e-07 3.536464e-07 4.715285e-08 0 7.072927e-08 5.894106e-08 0 4.833167e-07 0 3.064935e-07 3.536464e-08 1.768232e-07 4.715285e-08 2.357643e-08 1.414585e-07 5.316484e-06 2.357643e-08 8.841159e-07 5.894106e-08 5.068931e-07 1.178821e-07 3.536464e-07 0 1.178821e-08 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0.03968767 0.05102699 0.06236633 0.221117 0.02834833 0.085045 0.1133933 0.17009 0.1304023 0 0 0.01133933 0 0 0 0.01133933 0.02267867 0.102054 0 0.005669667 0 2.477644 0.833441 0.306162 0.476252 0 0 0.002267867 0.03968767 0.002267867 0.0005669666 0 0 0.3685283 0 0 0 0 0.1077237 0 0 0
0 0 0 0 0 0 0 0 0.005156302 0.01546891 0 0 0.03388427 0.1009162 0.007734453 0.02872797 0.007734453 0 0 0.04419687 0 0 0 0.004787995 0.01436398 0 0 0 0 0.002209844 0 0.02688643 0.02578151 0.002946459 0.01068091 0 0.01841537 0.01841537 0.0003683073 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01118543 0.001242825 0 0 0 0
0 0 0 0 0 0 0 0 0 0.004405214 0 0 0 0 0 0 0 0 0 0 0.0002936809 0 0.002936809 0 0.002936809 0.002936809 0.002936809 0.002936809 0.0002936809 0.5822224 0.04184953 0.04992576 0.008076225 0.0002936809 0.004405214 0.01248144 0.004405214 0.008076225 0 0 0 0.005139416 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.153048 0.006876764 0 0
0 0 0 0 0 0 0 0 0 0.005679257 0 0.00561472 0 0.0004302468 0.0005378084 8.604934e-05 8.604934e-05 0.0001720987 0 0.0009895676 0 0 0 0 0.004625153 0.002366357 0.000129074 8.604934e-05 0 0.09355715 0.0526622 0.001807036 0.008647959 0 8.604934e-05 0.0001505863 0.000322685 0.0001720987 0.0002151234 0.01157364 0 0.005184473 0.001462839 0.00202216 0.008970645 0.001376789 0.006346139 0 0 0
0 0 0 0 0 0 0 0 0 0.00940202 0 0.0002350505 0 0.000940202 0 0 0 0 0 0.0001175252 0 0 0.0001175252 0 0.003290707 0.01621848 0 0 0 0.3890086 0.1390324 0.01245768 0.01798136 0 0 0.002585555 0 0.0007051515 0 0.09237484 0 0.05182864 0.000470101 0.02091949 0.02221227 0.1557209 0.239634 0.0001175252 0 0
0 0 0 0 0 0 0 0 0 8.805065e-05 0 0.000264152 0 8.805065e-05 0.000264152 0 0 0 0 0 0 0 0 0 0.006163546 0.001849064 0 0 0 0.1280257 0.2356235 0.0001761013 0 0 0 0 0 0.0007924559 0 0.04948447 0 0.001232709 0.2740136 0.009509471 0.08179906 0.06771095 0.02333342 0 0 0
0 0 0 0 0 0 0 0 0 0.0005734076 0 0.001146815 0 0.0005734076 0.001146815 0 0 0 0 0.0005734076 0 0 0 0 0.01135347 0.02683548 0 0 0 0.2955343 0.1834904 0.002637675 0 0 0 0 0 0.0001146815 0 0.04013854 0 0.01559669 0.0974793 0.09736463 0.03830363 0.149086 0.1834904 0.001146815 0 0
0 0 0 0 0 0 0 0 0 0.04048849 0 0.01349616 0 0.09447315 0 0.0002699233 0 0.01349616 0 0.1079693 0 0 0 0.006748082 0.3657461 0.002699233 0 0 0 0.4588696 0 0.01349616 0 0 0 0.04048849 0.03009645 0.07422891 0 0 0 0.01430594 0.01349616 0 0.01349616 0.04453734 0.0008097699 0 0 0
0 0 0 0 0 0 0 0 0 0.01170313 0 0 0 0.01287344 0 0 0 0 0 0 0 0 0 0 0.1170313 0.02200188 0 0 0 1.755469 0.003510938 0.1165631 0 0 0 0.03510938 0.03510938 0.07279345 0 0 0 0.08192188 0.0002340625 0 0 0.008192188 0.0681122 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.04679543 0.04798313 0 0 0 0 0 0 0 0 0 0 0 0.08717727 0.4801876 0.0353935 0.439687 0.05059607 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.007798372 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0008664858 0 0.0181962 0.005198915 0.5718806 0.009531343 0.2521474 0.0008664858 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06391704 0 0 0 0 0 0 0 0 0 0 0 0 0 0.001620432 0 0 0.09209456 0.1186516 0.05068351 0.5731828 0 0 0 0
0 0 0 0 0 0 0 0 0 0.02056391 0 0.0008966821 0 0.0004782304 0 0.002749825 0 0 0 0.0298894 0 0 0.02492776 0 0.06695226 0.01661851 0 0 0.000597788 0.07478328 0.07478328 0.0004782304 0 0 0 0.004602968 0.03532927 0.0004782304 0 0.08966822 5.97788e-05 0.008129918 0.0298894 0.03353591 0.02152037 0.05858323 0.001374913 0.0007771245 0 0
0 0 0 0 0 0 0 0 0 0.0006562 0 0 0 8.2025e-05 0 0 0 0 0 0.00173893 0 0 0.000410125 0 0.00114835 0.0003281 0 0 0 0.00495431 0.0016405 3.281e-05 0.0022967 0 0.0006562 0.0016405 0.00082025 0.003281 0 0 0 0.00574175 0 0.00082025 0 0.00082025 0.00574175 0 0 0
0 0 0 0 0 0 0 0 0 0.004772 0 0 0 0 0 0 0 0 0 0.016702 0 0 0 0 0.002386 0 0 0 0 0.016702 0.009544 0 0 0 0 0.555938 0 0.028632 0 0 0 0.08351 0 0 0 0 1.283668 0.384146 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.002026578 0.0002026578 0.0002026578 0.0002026578 0 0 0 0 0 0 0 0 0 0 0 0 0.001215947 0.0002026578 0 0 0
0 0 0 0 0 0 0 0 0 0.01736132 0 0.02178232 0 0.000613132 0.003678792 0.001226264 0 0.001226264 0 0.007486664 0 0 0.04201568 0 0.04756614 0.03298005 0 0 0 0.03924045 0.05653722 0 0 0 0.006873532 0.01200448 0.002452528 0 0 0.003452901 3.227011e-05 0.006905803 0.001000373 0.001032643 0 0.003356091 0.01306939 0.0008390228 0 0
0 0 0 0 0 0 0 0 0 0.04643556 0 0.05825864 0 0.001713489 0.009766889 0.00325563 0 0.00325563 0 0.02004783 0 0 0.1124049 0 0.1273123 0.0882447 0 0 0 0.1050369 0.1511298 0 0 0 0.01833434 0.0322136 0.07110981 0 0 0.1009245 0 0.2009923 0.02947202 0.02998607 0 0.09749754 0.3814227 0.02484559 0 0
0 0 0 0 0 0 0 0 0 0.0008762082 0 0.0005154165 0 0 0 0 0 0.000154625 0 0.009535206 0 0 0.000154625 0 0.0002061666 0.0006700415 0 0 0 0.01577175 0.006184999 0 5.154165e-05 0 0.004123332 0.02798712 5.154165e-05 0.001958583 0 0.04509895 5.154165e-05 0.04040866 0.1643663 0.04499587 0.01396779 0.09612519 0.03834699 0.003762541 0 0
0 0 0 0 0 0 0 0 0 0.005705628 0 0.003169794 0 0 0 0 0 0.001267917 0 0.05895816 0 0 0.001267917 0 0.001267917 0.00443771 0 0 0 0.09762964 0.03867148 0 0 0 0.02535835 0.1737047 0.0006339587 0.0259923 0 0.06339587 0 0.5337932 2.170041 0.8495047 0.4615219 1.269185 0.506533 0.04944878 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.007884 0.031536 0 0 0 0.07884 0.03942 0.01971 0.007884 0 0.003942 0 0 0 0 0.03942 0 0.07884 0.3942 0.7884 0.1971 1.860624 0.1971 0.1971 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0147 0.000147 0.000294 0.000294 0 0.000294 0 0 0.000294 0 0 0 0.00147 0.00147 0.00147 0 0.00147 0.117747 0.00735 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0012168 0.024336 0 0 0 0.048672 0.006084 0.006084 0.0024336 0 0 0 0.0006084 0.0006084 0 0.0073008 0 0.0146016 0.0024336 0.0018252 0.0006084 0.0024336 0.0024336 0 0 0
0 0 0 0 0 0 0 0 0 6.94995e-05 0 0 0 3.474975e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01525198 0.005784254 0 0.003648724 0 0 0.006078047 0 0.000138999 0 0 0 0.0001895441 0 0 0 0 0.0003948835 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.193312 0.193312 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.97279 21.30298 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.04209933 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.85917 0.6701526 0 0 0.000687336 12.40401 2.52596 0 0.687336
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.1488804 0 0.007835808 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.507459 1.253729 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.34645 1.777288 1.340761 0 0
0 0 0 0 0 0 0 0 0 3.673133e-05 0 0 0 0 0 0 0 0 0 7.346265e-05 0 0 0.0004407759 0 0.0004407759 0.0004040446 0 3.673133e-05 3.673133e-05 0.002754849 0.003158894 0 0.004224102 0 0.004224102 0.004224102 0.00207532 0.00207532 0 0 0 0.0004407759 0 0.000110194 0 0.0006978951 0.01083574 0.0004775072 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.210228 0.4034093 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.04794171 0 0 0 0 0 0.2397086 0 0.1917669 0 0 0 0 0 0 0 0 2.157377 2.157377 0 0
0 0 0 0 0 0 0 0 0 0.0846932 0 0 0 0 0 0 0 0 0 0 0 0 0.0005081592 0 0.0005081592 0.0005081592 0 0 0 0.6080972 0.2134269 0.01693864 0.09824411 0 0 0.0846932 0 0.1693864 0 0 0 0 0 0.05589751 0 0.05589751 0.3048955 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.673244 2.673244 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0002784 0 0.005568 0.0843552 0.0843552 0.033408 0.033408 0.005568 0.002784 0 0.0286752
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.14565 0.14565 0.2913 0.72825 0.43695 0.5826 0.43695 0 0.14565
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0001399364 0 0.0002798728 0.004198091 0.004198091 0.001679237 0.001679237 0.0002798728 0.0001399364 0 0.001399364
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.2695643 0 0.2695643 2.156514 4.313029 1.527531 9.43475
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 748.7041
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 233.331
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 266.6193
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 97.5
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19.80945 190.1707 146.5899 39.6189
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 850.2666 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
3.42924e-06 0.6171009 2.463139e-06 1.157301 0.07552549 0.0025284 0.1498636 0.0325392 0.05007185 0.2380649 0.3492465 0.246516 0.4554724 0.4866565 0.2764926 0.1851866 0.1873525 0.1421396 0.006657664 0.5056064 0.00084103 0.06696411 0.3548135 0.1071614 1.313005 0.8336864 0.03061366 0.02518649 0.006555066 8.184005 3.524776 0.7970927 0.7363228 0.00757444 0.3355124 1.057455 0.3514416 1.109396 0.05886494 0.6191112 0.002858881 3.69705 212.852 63.90867 75.55833 32.76805 161.9631 367.6832 884.8559 0
];

