food_web_filename = 'Mhlanga_Sep_2002';
crc = 'crc_bc02c773';
n = 46;
l = 43;
b0 = [0.9670837 1.788513 0.8624786 2.955855e-07 0.00790624 4.923839e-05 0.001805304 9.482482e-06 0.0002614428 3.954264e-05 0.0004066845 0.0008775684 0.0001179142 0.06795349 0.1924161 0.01424716 0.005205167 0.0004722188 0.07630205 0.08794952 0.008432876 0.0570118 0.7347405 0.6772055 0.02028755 0.2221522 0.002615109 24.61887 6.411409 0.01451472 0.009559806 0.06549914 0.04030975 0.02136105 0.2729114 0.68095 0.4539769 0.5547164 0.4470768 0.007402939 2.291469 0.204617 0.1338769 81.1918 1450.554 116.8958].';
p = [2295.1190556 174.83962482 32.255688654 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 906.3615456 0 1842.7260852].';
q = [0 0 9.935876412 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.067358293 0.5119323426 0.0011816319816 0.0007782979806 0.005349682548 0.003251901051 0.0017459824158 0.02253887685 0.05608539216 0.03739711878 0.04576017474 0.0365062005 0.0006036163182 0.18861502212 0.016898445522 0.01615245219 913.874535 0 1328.738733].';
r = [2077.0629516 68.68215522 223.9394094 0.000231071652 1.0074618792 0.006233881122 0.12341976948 0.0012015999954 0.03306424212 5.323611258e-06 0.0515327274 0.1110530988 0.01321691994 69.46539012 194.5930392 11.348943984 5.27031036 0.3780249354 46.23782562 28.42565355 3.834518562 19.286904546 9.300076758 8.560488888 0.254672838 2.791442241 0.003809872332 183.1150944 22.49683758 0.0519421581 0.03425321718 0.2356177761 0.14284048716 0.07683090408 0.9884903112 2.465823339 1.6441838028 2.0140913772 1.6054875774 0.026590172616 8.302286286 0.741299769 1.0167390198 0 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1010.6171082
7.506331812e-05 0 0.0004663231146 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.13947763788 0 0 2.859337215e-05 0 0 0.0614352501 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.5418769642 0 0
0.0010520392176 0 0 2.908744776e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.018511770564 0 0
0.022213447956 0 0 2.9262759372e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3840383232 0 0
0.00020166265224 0 0 2.7268311168e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.003538694376 0 0
0.005551551432 0 0 2.719191699e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.09851253426 0 0
1.1182361778e-06 0 0 3.5916879726e-06 8.777244168e-06 7.733833506e-06 0 7.64705466e-06 7.750123038e-06 0 7.37455761e-06 7.587264258e-06 6.193135368e-06 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0410838406e-05 0 0
0.00606129174 0 0 2.7214034652e-05 0 0 0.04088723688 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.1072987965 0 0
0.014480739882 0 0 2.7198282258e-05 0 0 0.06185023488 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.25677264564 0 0
0.0008387311212 0 0 2.7191690442e-05 0 0.0026526636864 0.005670233352 0.0005069417094 0.005666729166 0 0.005393318616 0.005548167828 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.014873468148 0 0
0 9.058635054 0 0 0 0 0 0 0 0 0 0 0 8.33018823 0 0 0 0 0 4.599322686 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 82.37661642 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 47.26131102 2.3555478114 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 241.31037924 0
0 1.655561628 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 15.236916912 0
0.33243112728 0.18924058692 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.615671824 1.752927939 0
0 0.0563987592 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5235827562 0
0 0 0 0 0 0 0 0 0 0 0 0 0 8.294130558 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 60.0438069 0
0 2.4215445324 0 0 0 0 0 0 0 0 0 0 0 8.34074829 0 0 0 0 0 4.597850124 0 4.709764836 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 21.886425918 0
0 0.554540175 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.149483542 0
0 2.841447231 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 25.971586956 0
0 11.326976136 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 125.30994714 0
0 11.00012949 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 120.04186572 0
0 0.456121512 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.551265278 0
0 4.957988364 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 49.32720954 0
0 0.0008878467276 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00875441763 0
0 27.06263847 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 273.09670794 0
0 1.2318744816 0 0 0 0 0 0 0 0 0.019994662422 0.04525782282 0.005535305748 0 0 0 0 0.07264583172 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 43.77882978 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0.005949472032 0.010831358538 0.00595041111 0.01223534592 0 0.009924907104 0.005881635522 0.00594379737 0.00595134288 0.007900729956 0.007595109396 0.005949884934 0.00596824263 0.0006803101242 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.010914256836 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0.003696433902 0.0067292064 0.003697014888 0.007601405238 0.0035514910494 0.006166092114 0.00365535198 0.003692253726 0.003696938154 0.00490774221 0.00471791691 0.003696031962 0.00370743975 0.000680394897 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.006779507364 0
0 0.003944332224 0 0 0 0 0 0 0 0 0 0 0 0.023753824542 0.0432502056 0.023757569892 0.04885719552 0.022822054542 0.0396301878 0.023489713422 0.023736307266 0.023766452766 0.031553254656 0.030332453256 0.023760624636 0.023833960416 0.0006801292512 0.028503626256 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.04359072186 0
0 0 0 0 0.06514632558 0.0026337304854 0 0.0005032825938 0.01824317964 3.4856233146e-05 0.019996372494 0.04526158644 0.005535872118 0 0 0.0256810428 0.052805781 0 0 0 0 0 0 0 0 0 0.0006803115858 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.047114676 0
0 0.0012849806214 0 0 0 0 0 0 0 0 0 0 0 0.007728462126 0.014070099708 0.007729682562 0.015893898804 0.00742538475 0.012892605516 0.007642556586 0.007723328256 0.007733131938 0.010266146856 0.009866826774 0.007729510824 0.007753364136 0.0006803115858 0.009272046924 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.014178733128 0
0 0.10592858304 0 0 0 0 0 0 0 0 0 0 0 0 0 0.637562709 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.1698167726 0
0 0.3958728984 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.371674832 0
0 0.26492110218 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.914846956 0
0 0.32339931624 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.5712160722 0
0 0.031008435948 0 0 0 0 0 0 0 0 0 0 0 0.18690842142 0.34034426622 0.18693798228 0.3844713222 0 0.31185494172 0.17776121706 0.17964072504 0.1798688808 0.23881122594 0.22957044264 0.17982477702 0.18037992924 0.0006800031882 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3299295258 0
0 0.0004969381536 0 0 0 0 0 0 0 0 0 0 0 0.0029871512118 0.005438007036 0.0029876244048 0.006142852674 0 0.00498294153 0.0029539531602 0.0029851674552 0.0029889584802 0.003965602158 0.003812214546 0.0029864989728 0.0029957139954 0.0006803897814 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.005478047568 0
0 1.333291617 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 14.712034428 0
0 0.1194921945 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.3146558516 0
0 0 0 0 0.7766032554 0.0026059516812 0 0.0004986526104 0.018073754622 0 0.0198055935 0.044841888 0.005469895494 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.8168042438 0 0
159.48124164 0 404.5405518 4.268756268e-05 0.4483881864 0.0027315725472 0.05651663724 0.0005248031922 0.014521091004 1.9005703668e-05 0.01877417892 0.04058446644 0.005707339722 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11.07361143 0.02428227333 0.015975996876 0.10917103302 0.06878136132 0.0356548185 0.4517717904 1.1231160948 0.7498048194 0.9181376568 0.748797777 0.012346770996 3.778718328 0.33806475486 0.82595016 0 36.04413393 0
0 17.490945276 0 0 0 0 0 0 0 0 0 0 0 9.703281042 48.65213304 2.293685226 2.0919555594 0.095511906 21.714673302 4.112112942 1.6457831586 4.592361816 127.03945842 122.1956064 4.528763946 51.2690913 0.0010705427082 114.9391278 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 558.7784496 0 0
58.05266922 13.272910182 404.4565098 4.102751394e-05 0.445210668 0.0027273697164 0.05650176546 0.000524697957 0.014514543036 1.8998421246e-05 0.018770258178 0.04057624494 0.005705914662 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11.071762506 0.024270774192 0.015967694988 0.10912407912 0.06876323748 0.035639387658 0.4505067756 1.122525243 0.7483819518 0.9166263624 0.747200979 0.012341487312 3.775703778 0.33788475882 0.8258617332 0 0 0
];

