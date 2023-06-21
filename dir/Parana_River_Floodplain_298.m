food_web_filename = 'Parana_River_Floodplain_298';
crc = 'crc_cfcde9e2';
n = 40;
l = 39;
b0 = [4.8 0.6000739 0.08442883 0.001 0.004088625 0.00259016035 0.00029065102 0.050501518 0.00142387743 0.00136458932 0.005 0.000436241215 0.04 0.002532565 0.006 0.0208380334 0.34 0.0314646475 0.06 0.4460186 0.016 0.0145169077 0.006 0.001874394 0.09 0.00691664871 0.01 0.04770644 0.006 0.00213029166 0.0296022743 0.0249465611 0.0146160843 0.0254259948 0.00281198532 0.4135749 0.06337104 35 7.4 17.0240441506934].';
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 19.011312 420 177.6 0].';
q = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 264.11759389048].';
r = [103.679993 105.012924 12.2421808 0.01754 0.04040379 0.04281535 0.008568392 1.66857016 0.0337743722 0.04589114 0.06573 0.008352274 0.59384 0.0283140764 0.11442 0.265809923 7.191 0.56825155 0.707999945 6.155057 0.212640017 0.0548158437 0.024120003 0.008491005 0.273420036 0.03582824 0.06599999 0.179757863 0.017975999 0.004192414 0.112015009 0.481468618 0.439944148 0.4500401 0.02249588 8.850503 3.168552 70 29.6 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 96 96
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 45.00554 0 105.0129
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10.5536 0 0 10.5536
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.005 0 0 0.02
0.0151688 0.00606752 0.00606752 0 0 0 0 0 0 0 0 0 0 0.00606752 0 0.00303376 0 0 0 0.00606752 0 0 0 0 0 0 0 0.00303376 0 0 0 0 0 0.00303376 0 0 0.00303376 0.00303376 0 0.00606752
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0194262 0.04532781 0 0
0.001220734 0.002441469 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.001220734 0 0 0.007324406
0.3560357 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3560357 0 0 1.6615
0.004812706 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.01684447 0 0 0.02646988
0.006550029 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0229251 0 0 0.03602516
0.0337225 0.01927 0.009635 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.009635 0.0240875
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.003431473 0.008006771 0 0
0 0 0.58436 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.25044 0 0 0
0 0.006458041 0.01291608 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.02367948
0.0348 0.0348 0.0348 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0174 0 0.0522
0.01976487 0.0790595 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03952975 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03952975 0.01976487 0.1185892 0 0.0790595
0.51 2.04 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.02 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.02 0.51 3.06 0 2.04
0 0.4946243 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3297495 0 0
0.0963 0.0963 0.0963 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0963 0.4815 0 0.0963
0 2.809917 0.9366391 0 0 0 0 0 0 0 0 0 0 0 0 0 0.9366391 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.9366391 0.9366391 1.873278 0 0.9366391
0 0.224 0.064 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.032 0 0 0 0 0 0 0 0
0 0.04686784 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.008521425 0 0 0 0 0 0 0 0 0 0 0 0 0 0.004260712 0 0 0 0.02556427 0 0
0 0 0 0 0 0.00468 0 0 0.00468 0 0.00936 0 0 0 0 0 0 0 0 0.00936 0.00468 0 0 0 0 0 0 0 0 0 0 0 0.00468 0.00468 0 0.00468 0 0 0 0
0 0 0 0 0 0 0 0.003373909 0 0 0 0 0 0 0 0 0.003373909 0 0.003373909 0 0 0 0.001686955 0.003373909 0 0 0 0 0 0 0 0 0 0.001686955 0 0 0 0 0 0
0 0 0 0 0 0 0 0.14067 0 0 0 0 0 0 0 0.04689 0 0.070335 0 0.070335 0 0 0 0 0 0 0 0.04689 0 0 0.04689 0 0 0.04689 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.02074995 0 0.006916649 0 0 0.0138333 0.006916649 0 0 0 0 0 0 0.006916649 0.006916649 0 0 0 0 0 0.006916649 0 0 0 0 0
0 0 0 0 0 0 0 0.02 0 0 0 0 0 0 0 0 0 0.02 0 0.02 0 0.01 0.01 0 0 0 0.01 0 0 0 0 0.01 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0.03439634 0 0 0 0 0.03439634 0 0 0 0 0 0.03439634 0.03439634 0 0 0 0 0.03439634 0.01719817 0 0.03439634 0 0 0 0.03439634 0.03439634 0.03439634 0 0.01719817 0 0 0 0
0 0 0 0 0.008055 0 0 0 0 0.006444 0 0 0 0 0 0 0 0 0 0.004833 0 0.003222 0 0 0 0 0 0.003222 0 0 0 0 0.003222 0.003222 0 0 0 0 0 0
0 0 0 0.001580676 0 0 0.001185507 0.001185507 0 0 0 0.0007903382 0 0 0 0 0 0 0 0.0007903382 0 0 0 0.001580676 0 0 0 0 0 0 0 0 0 0 0 0.0007903382 0 0 0 0
0 0 0 0 0 0 0 0.01992233 0 0 0 0 0.01992233 0 0 0 0 0 0.01992233 0.01992233 0.01992233 0 0 0 0.03984466 0 0 0 0 0 0 0.01992233 0.01992233 0 0 0.01992233 0 0 0 0
0.2532076 0.07234503 0.07234503 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03617251 0.07234503 0 0.2170351
0.06284916 0.06284916 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06284916 0 0 0.4399441
0 0.4119011 0.06865019 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06865019 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.1373004 0 0
0 0 0 0 0 0.004217978 0 0.008435956 0 0 0 0 0.004217978 0 0 0 0 0 0 0.008435956 0 0 0 0 0 0.002108989 0 0 0 0.002108989 0 0 0 0 0.004217978 0.008435956 0 0 0 0
2.047196 4.094392 1.364797 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.364797 2.047196 1.364797 1.364797
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
84.87837 34.50425 5.614517 0.005879324 0.0122164 0.01304068 0.002453443 0.4770171 0.009672684 0.01316515 0.02126 0.002295632 0.1824233 0.008672009 0.03883005 0.07956378 2.06207 0.1657872 0.1973074 1.886659 0.07584102 0.01717641 0.01099305 0.003423956 0.121239 0.01402835 0.024 0.06974682 0.007327351 0.001601979 0.0403183 0.145663 0.1263268 0.138292 0.008548435 2.750273 1.584276 296.7752 50.62556 0
];

