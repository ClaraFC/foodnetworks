food_web_filename = 'Mhlanga_Dec_2002';
crc = 'crc_f9f9e69f';
n = 47;
l = 44;
b0 = [0.3345228 7.584262 0.05847582 3.96409e-08 0.009902149 0.003431783 2.658357e-05 2.34805e-05 0.001620243 4.710133e-05 0.00168984 1.065213e-06 0.003034111 0.0364007 0.0009014574 0.0006067538 7.087313e-05 0.009231564 1.259514 2.819751 0.01631808 0.05509019 56.1385 0.02549298 11.86394 0.09727074 0.04837764 0.00476112 1.078519 0.1093593 0.004881487 0.0004698556 0.04342848 0.1464869 0.8288417 0.05440025 0.3742383 1.168618 0.293534 0.002280609 0.001277488 0.003408012 0.005171534 0.0172164 2.053942 1185.292 4.17336].';
p = [91.22360814 3931.532262 88.25553702 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2768.3975592 0 5698.94283].';
q = [0 0 76.83072138 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0485434358 0.008064849366 0.003874036572 0.0003773551572 0.0859079151 0.008723534022 0.0003937780602 3.397352175e-05 0.0034619205474 0.011752751178 0.06611898384 0.004337020296 0.02981786409 0.09289081872 0.02352200382 0.0002777776281 0.0001545667578 0.0004131552222 0.0006638207184 0.002221231887 2844.8575092 0 5783.562162].';
r = [35.536103694 2896.149438 22.889152944 2.7158581548e-05 1.3138044696 0.4497218964 0.001601643204 0.0025653827808 0.00017330549292 0.005575148964 0.17969937174 7.40380788e-05 3.092301639 28.5211143 1.0015040322 0.7937043408 0.06445096938 2.804523561 17.814349854 38.92686588 0.21868883064 0.7683675048 713.1025818 0.03707381286 92.87411994 0.3544396443 0.17015369868 0.016541592228 3.777709824 0.3829841442 0.017334791586 0.0014977566954 0.15260715414 0.5151438432 2.91238416 0.19019497518 1.3110869898 4.089732192 1.034968095 0.017467238124 0.009702729288 0.02599121259 0.0419886621 0.13933184328 0 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 106.38763506
6.588421434e-05 0 4.43484153e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1.2816766746 0 0 3.80458134e-06 0 0 0.0005107974102 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.5425317152 0 0
0.4502382066 0 0 3.884154498e-06 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5330020374 0 0
0.0022573582542 0 0 4.612977684e-06 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0026571080466 0 0
0.003603297096 0 0 4.705613892e-06 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.004233009186 0 0
0.00028194234768 0 0 3.277555785e-06 0.00030582361278 0.00030373378056 0 0.0002630620566 0 0.00027656464626 0.00027778588614 1.534360275e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00032719118166 0 0
0.006555798522 0 0 2.896467336e-06 0 0 0.0005795726328 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00778607109 0 0
0.14364469602 0 0 2.8525831614e-06 0 0.15736488138 0.0005711907222 0.0010826717958 0 0.0011292745464 0 1.3431223386e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.17101776006 0 0
6.288428034e-05 0 0 2.9335072032e-06 0 0 5.989348134e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.265566098e-05 0 0
0 1.0410826986 0 0 0 0 0 0 0 0 0 0 0.2056078395 0 0 0 0 0.1868977152 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.2369097852 0
0 0 0 0 0 0 0 0 0 0 0 0 0 5.695533648 0.0552152286 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 36.258484878 0
0 0.3718672146 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.1519319042 0
0.14328181728 0.2149204239 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.17275191192 0.6695472672 0
0 0 0 0 0 0 0 0 0 0 0 0 0.03945008214 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06093165582 0
0 0.9291731022 0 0 0 0 0 0 0 0 0 0 0.2053533384 0 0 0 0 0.18666846324 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.8938550914 0
0 64.27181376 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 246.3442758 0
0 111.55252752 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 432.7187382 0
0 0.8963163342 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.109543038 0
0 2.8710402462 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10.729883304 0
0 540.972873 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2167.258653 0
0 0.02225830446 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.06992539218 0
0 39.45925368 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 134.33075964 0
0 0 0 0 0.25543922796 0 0 0 0 0.0010174073364 0.06141135294 1.0726587396e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.27956428794 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.018589534992 0.028803246948 0.018234991026 0.012439557018 0.0029827141596 0.018753511896 0.018975667788 0.023575414338 0.018963134568 0.018996111918 0.08615741022 0.004531230396 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.028817135802 0
0 0 0 0 0 0 0 0 0 0 0.0004437881658 0 0.0013555287648 0 0.0013301813322 0 0.0016449091218 0.001367714124 0.0013837950126 0.0017167438386 0.0013828877244 0.0013852748826 0.0062466957 0.0013638989826 0.0013830960024 0 0 0 0 0 0 0 0 0.0013776983136 0.001981794402 0.0018124618302 0.0018854724042 0 0 0 0 0 0 0.0004960348848 0 0.0020893564692 0
0 0.21452729004 0 0 0.15696484146 0 0 0 0 0.0010189113228 0.06150021822 1.0743110784e-05 0.20372208318 0.6751831968 0.05548244562 0.06323144688 0.0029826019818 0.18519443964 0.4448058048 0.552632787 0.4445123886 0.4452852096 2.0196732276 0.004531025772 0.4445796222 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.6755018256 0
0 0.019051689258 0 0 0 0 0 0 0 0 0 0 0.03869118288 0.05996546514 0.03796206372 0.02589584319 0.0029824857846 0.03904163802 0.03950419788 0.04908096648 0.03947810832 0.03954676698 0.1793781486 0.004530810186 0.039479643 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.05998760838 0
0 0 0 0 0.0019720835316 0.001969277625 0 0.00097784694 0.0005680347624 0.0010246006008 0.0017633458584 1.0805630724e-05 0 0 0.005349247722 0.0036521678844 0 0 0 0 0 0 0 0.004539082842 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.008448018768 0
0 0 0 0 0.0005665128714 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0014142463524 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0020949010488 0
0 0.007426945008 0 0 0 0 0 0 0 0 0 0 0.01507965606 0.023369383548 0.014795564868 0.010093842486 0.0029830185378 0.015211817586 0.015391996326 0.019122595128 0.015381830898 0.015408578178 0.06987942486 0.004531789458 0.015384165804 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.023373863352 0
0 0.20437158168 0 0 0 0 0 0 0 0 0 0 0 0 0.05546804886 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.6436897362 0
0 0.18226926648 0 0 0 0 0 0 0 0 0 0 0.2037357126 0.5736038238 0.05548193406 0.0632377683 0.0029827660464 0.18520675362 0.377894853 0.469495152 0.3776452848 0.3783019086 1.7157711438 0.0045313254 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5738804316 0
0 0.00995724135 0 0 0 0 0 0 0 0 0 0 0.020217019284 0.031330824714 0.019836144594 0.013529749842 0.0029830349808 0.020395718154 0.02063729871 0.025639182576 0.02062366929 0.020659529646 0.09369254286 0.00453181869 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.031339166796 0
0 0.0731560032 0 0 0 0 0 0 0 0 0 0 0.1486238922 0.230373738 0.05545522332 0.06319315296 0.0029763002934 0.14997072006 0.15174809874 0.18854837316 0.15164330202 0.15190712082 0.6892160184 0.004529198772 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.23047795008 0
0 0.26899739496 0 0 0 0 0 0 0 0 0 0 0.20375496918 0.8464023288 0.05548851126 0.06324876684 0.0029777425272 0.18522809298 0.5576380362 0.6927925536 0.557269713 0.5582387538 2.4050492802 0.004531847922 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.804466467 0
0 0.443471364 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.354561551 0
0 0 0 0 0.0005121179658 0 0 0 0 0.0004539908646 0.0004565888586 1.042036758e-05 0 0 0 0 0 0.0014585083506 0.0014759463348 0.001837000998 0.001474962678 0.0014775511716 0.006749350902 0.0014543712918 0 0.0005284243062 0.0020554217712 0.0017420273262 0.0018895824234 0.0019162500462 0.000942452469 0.0009662725296 0.0019087487496 0.0014742019152 0.0021314267982 0.0019472012532 0.002026632636 0.0021464885862 0.0014682060666 0 0.0021831984972 0.0020321534646 0.000507991869 0.0005134560606 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0019067507424 0.0019295140662 0.0023953234284 0.001923807249 0.0019271787948 0.00879313176 0.0018969890814 0.001924101396 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0029311898364 0
0 0 0 0 0 0 0 0 0 0 0 0 0.004208337378 0 0.00412834401 0 0.0029559518982 0.004246795728 0.00429312114 0.005342783796 0.004290260058 0.004297787298 0.019624004316 0.004230396576 0.004290917778 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.006538946274 0
0 0 0 0 0.034391133756 0 0 0 0 0.0009900023364 0.030663629892 1.0425435678e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.03771892656 0 0
0 0 0 0 0.09130132872 0.09124118388 0.0004976525106 0.000940214394 0.0005421534804 0.0009875782728 0.05961946788 1.0398794364e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.10015347258 0 0
26.742116898 0 47.32554834 5.118398964e-06 0.4838926428 0.14091640038 0.0005472986778 0.0010024065702 0.0003840697476 0.0012214423116 0.0393678306 1.6233969276e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0.11744770842 0.0618271416 0.005988314052 1.3918758336 0.14037268518 0.00579292371 0.000787659894 0.05468908914 0.18671158044 1.0923680502 0.0684922203 0.4731714414 1.5120606438 0.3703106106 0.013986901782 0.006790787262 0.019992992544 0.030219826014 0.10107665568 0 116.83544418 0
0 161.1790362 0 0 0 0 0 0 0 0 0 0 0.26980777908 5.323552794 0.08806695408 0.08827541478 0.007479259326 0.22897732536 291.16613988 503.3121912 2.1525816312 11.1951252 1987.8297138 0.005376100968 79.3588509 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 117.9123876 0 0
26.913697776 110.1776004 47.5977348 4.639845546e-06 0.4855720212 0.14172689412 0.000551030508 0.0010094259042 0.000387161397 0.001229419359 0.03962339136 1.5800093316e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0.11696220144 0.0619093566 0.005998044654 1.3939582482 0.14057972082 0.005810561568 0.0007899973578 0.05476757706 0.18706923396 1.089052776 0.0685890513 0.473830623 1.5092532756 0.367764138 0.014009037714 0.006796699434 0.020018117448 0.030393833148 0.10165420692 0 0 0
];
