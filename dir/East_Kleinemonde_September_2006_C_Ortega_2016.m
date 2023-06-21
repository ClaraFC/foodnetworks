food_web_filename = 'East_Kleinemonde_September_2006_C_Ortega_2016';
crc = 'crc_9eab0606';
n = 45;
l = 43;
b0 = [0.137 0.281 0.0322 0.39 0.00627 0.0276 0.0253 0.0685 0.00243 0.0263 0.107 0.000276 0.000914 0.0157 0.00156 0.0281 0.0214 0.00873 0.00643 0.00336 0.158 0.0249 0.00873 0.862 0.0198 0.0202 1.11 0.174 0.00106 1.35 0.00228 0.163 0.925 0.00568 0.0805 0.00143 0.615 0.109 0.0104 0.00958 0.0279 0.00181 0.00189 0.622 71.7].';
p = [15.96798 32.92254 3.76362 9.64656 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 36.13806 0].';
q = [0.642794495173492 1.48243273595806 0.0742098580637633 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 14.6826585342775 0].';
r = [11.9934709105139 27.7120282709044 1.3921199370898 3.40193984815387 0.942933831656658 0.84419390371271 0.721696184102167 0.767922173679136 0.181636671174176 1.19561622098859 2.95804582936484 0.0266455050564275 0.0723352824562379 0.855353950580722 0.11267603949725 1.24560616310665 0.998638458083041 0.493117151259768 0.338470618130287 0.206667195587964 1.12792799535505 1.92085112574876 0.493117184243083 3.63667766529034 0.234873180042337 0.111812285065988 2.803863310987 1.14439960172385 0.00435356790851905 1.81377723797302 0.0159371478441833 0.00322523375265671 0.119466190293103 0.00492737398981171 0.25780580147037 0.0115252921700565 0.115639364951232 7.26736829564321 0.597859540778414 0.412034638063262 2.74186660025859 0.163174795124305 0.0930668027510795 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.277581977563266 0.310042932329409 0.135996598578694 0.134506343512937 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.434852256241667 0
0 0 0.0841483521028296 0.0832397311932816 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.772702958811214
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.808679027783873 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.974215222806961 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.703437473332108 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00647567643528765 0 0 0 0 0 0 0.474090182701454 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.16256777515814 6.21013399725758
0 0 0 0.260973731423361 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.842297857678187 2.42008812218025
0 0 0 0.00197137115797961 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00636930657671618 0.0183050918530112
0 0 0 0.00599267599459973 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0193673693956373 0.0528903113244951
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00647567643528765 0 0 0 0 0 0 0.474090182701454 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.24718550466623 3.58316324750536
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.562749978665687 0
0 0 0 0.137210355851963 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.27382677930798
0 0 0 0.109183633365024 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.01102004922981
0 0 0 0.0534225140494195 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.493082446280193
0 0 0 0.0277599753237647 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0897111167288481 0.254984666067504
0 0 0 0.0167283357946814 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.05404592327575 0.155308177482916
0.170145035220454 0.190447114119137 0.0832348875771802 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.266063804848368 0.764296067584791
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.47500363993681 7.10930584147498
0 0 0 0.0534225140494195 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.493082446280193
0 0 0 0.857975807031717 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 28.7294826770755 0
0 0 0 0.119597734941006 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.385520277592402 0
0 0 0 0.0787452241973318 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.253529828050814 0
0 0 0 0.857975807031717 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 29.4091677162432 0
0 0 0.997137876198967 0.857975807031717 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.91858098236373 0
0 0 0 0.000995368866420102 0.000546024474610483 0.000372760944496521 0.000322296726267398 0.000381769111082158 0.000268509560162952 0.00195078073109634 0.00114659455909026 8.42744954877546e-10 2.57175870319399e-05 0.00129509360978953 0.00021165563983357 0.000548607644020857 0.00044119934806231 0.000229707055486759 0.000147434453464807 9.04146139495355e-05 0.000550993008420933 0.00192158218231422 0.000229670516575188 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0.857975807031717 0.0283245627726322 0.00492263717879229 0.00454942657826429 0.0108831596355383 0.00982343139153983 0.840696211639795 0.0608568677362221 8.42744954877546e-10 2.57175870319399e-05 0.32759325235489 0.00640993352731519 0.00905696030568743 0.00569867260670374 0.0014519357485573 0.000711223218746937 0.000235136469745882 0.0224745225119175 0.474090182701454 0.00145193584567348 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0.0137466140802952 0.00737900545005734 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.0283245627726322 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0.857610399958876 0.0283245627726322 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.00152477517274092 0.00104044157743294 0.00090001908933855 0.00106639620597973 0.000749634853679425 0.0054462493710578 0.00320183880702286 8.42744954877546e-10 2.57175870319399e-05 0.00361661196027592 0 0.00153215405978377 0.00123221458352779 0.00064152750219488 0.000411895461216751 0.000235136469745882 0.0015389804717964 0.00536595519048369 0.000641527545104975 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.0283245627726322 0.00492263717879229 0.00454942657826429 0.0107662542618098 0.00982343139153983 0.295230578322375 0.0604913610230917 8.42744954877546e-10 2.57175870319399e-05 0.196036890468565 0 0.00905696030568743 0.00556708042335458 0.0014519357485573 0.000711223218746937 0.000235136469745882 0.0224745225119175 0.290266009316287 0.00145193584567348 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.000701718200302629 0.000478741997343571 0.000414016089411521 0.000490637240366831 0.000345007385874683 0.00250673679105466 0.0014733575606289 4.21372477438773e-10 2.57175870319399e-05 0.00166428484486924 0 0.000705039403941528 0.000566942989929281 0.000295139205283274 0.000189464602568557 0.000116142943868886 0.000708106399416292 0.0024691435498445 0.000295139225024363 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.0283245627726322 0.00492263717879229 0.00454942657826429 0.0107662542618098 0.00982343139153983 0.840696211639795 0.0604913610230917 4.21372477438773e-10 2.57175870319399e-05 0.32759325235489 0 0.00905696030568743 0.00556708042335458 0.0014519357485573 0.000711223218746937 0.000235136469745882 0.0224745225119175 0.474090182701454 0.00145193584567348 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.89135106353282 0.00922977057851499 0.0044761453988835 4.00865118226604 0.977051256388754 3.91894158208133e-05 0.217507143502344 5.55114138392902e-05 0.0101464742717305 0.210545460140414 0.000236068332810624 0.112575443487062 2.52095217957401e-05 0.462411358269929 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.305644616398611 0.00922977057851499 0.00443229744395566 0.314991551423275 0.0768048519420031 3.91894158208133e-05 0.0170966313669965 5.55114138392902e-05 0.000808135973155482 0.0165499034259253 0.000236068332810624 0.0303025917965173 2.52095217957401e-05 0.0363391042924765 0 0 0 0 0 0 0 0
0 0 0 0 0.0283245627726322 0.00492263717879229 0.00454942657826429 0.0107662542618098 0.00982343139153983 0.149680477681085 0.0604913610230917 4.21372477438773e-10 2.57175870319399e-05 0.0989871153065211 0.00640993352731519 0.00905696030568743 0.00556708042335458 0.0014519357485573 0.000711223218746937 0.000235136469745882 0.0224745225119175 0.146905816829387 0.00145193584567348 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.4319441143648 0.00922977057851499 0.00443229744395566 1.47592909071764 0.359872020350518 3.91894158208133e-05 0.0801130610850925 5.55114138392902e-05 0.00378665250198589 0.0775616886826479 0.000228470899231077 0.112575443487062 2.52095217957401e-05 0.170062186738135 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.083161643385922 0.00289572176701233 0.00198667775785548 0.0856908570867262 0.0208966198028694 3.91894158208133e-05 0.00465042985687746 5.55114138392902e-05 0.000215247071276012 0.00440600077961721 0.000228470899231077 0.00806875907407266 2.52095217957401e-05 0.00967922669364385 0 0 0 0 0 0 0 0
0 0 0 0 0.00672479941956686 0.00458642142493269 0.00396841547308837 0.00470178799964343 0.00330650344086378 0.0240219804473771 0.0141195243282298 4.21372477438773e-10 2.57175870319399e-05 0.015948330283691 0.00260597582435293 0.00675434327881775 0.00543183290157906 0.0014519357485573 0.000667000192298643 0.000235136469745882 0.00678510621112515 0.0236643010239723 0.00145193584567348 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
2.88398758152882 3.22758894668906 0.996772490388707 0 0.163222577216226 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20.2386793811345 0.23965979898844 0.205135349137404 21.578017530794 4.19467031538628 0 0 0.00496642603078311 0.0101428192017479 0.457405719409681 0 0 0 0 2.63693298087768 0.214695892547261 0 0.98398810694248 0.0588247694021789 0 0 0
0 0 0 0.857610399958876 0 0.0697282237352316 0.0631803899905137 0.15647053614975 0.477837421351212 5.49742218394022 0.303041615856467 2.59474812646599e-07 0.00570933356223565 3.48282582954407 0.421485087779016 0.11966298644399 0.0914931208119272 0.0449617565640856 0.0297244524052869 0.0177978645890961 0.246777637856431 5.29650481676489 0.0449617595714597 0 0 0 0 0 0.00716585590311609 0.534087073900617 0 0 0 0.0233146242972336 0.420057624951723 0.00181920618088372 1.00810059148894 0 0 0.149800891019276 0 0 0.0333842455708346 5.20726441297794 0
];
