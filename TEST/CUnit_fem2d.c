#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "mkl.h"


//#define NDEBUG
//#define MATLIB_NTRACE_DATA

#include "fem2d.h"

/* CUnit modules */
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/TestDB.h>

/*============================================================================*/

static const matlib_real TOL = 1e-9;



/* Nodes */ 
#define SIZE_p 178
static const matlib_real p[SIZE_p] =
{
    -1.002597122608271e+00,    -3.093844724217988e-01 ,
    -1.001866477520648e+00,     7.922530463708422e-02 ,
    -9.984656279523078e-01,     5.237956526415085e-01 ,
    -9.848825573603598e-01,    -6.794737361577037e-01 ,
    -9.268331503612441e-01,     8.482382054115748e-01 ,
    -8.865693774220601e-01,    -8.966087502396579e-01 ,
    -7.992094181033553e-01,     2.769047151900140e-01 ,
    -7.688818284902823e-01,     9.632767089259799e-01 ,
    -7.590934064821773e-01,    -1.886569430247549e-01 ,
    -7.550980256978949e-01,    -4.600355102615336e-01 ,
    -7.468550817945419e-01,     7.059154350108767e-01 ,
    -7.329734295049922e-01,    -7.252725504867040e-01 ,
    -7.079354830257115e-01,     4.517263970311526e-02 ,
    -7.034503864537140e-01,     4.663234775057328e-01 ,
    -6.480564787465959e-01,    -9.897232463432050e-01 ,
    -5.946231948951380e-01,    -5.758793851907518e-01 ,
    -5.830002370108337e-01,     2.440401645415483e-01 ,
    -5.697800755496611e-01,    -1.380158074381450e-01 ,
    -5.640851908298194e-01,     5.859483497519962e-01 ,
    -5.489431992410468e-01,    -3.476148975392481e-01 ,
    -5.225366421750383e-01,     7.724559133344134e-01 ,
    -4.874108260327063e-01,    -7.627958847855207e-01 ,
    -4.852214751773601e-01,     4.915159352862503e-02 ,
    -4.690937877875496e-01,     4.107291010144715e-01 ,
    -4.096450489633402e-01,    -5.290934923931648e-01 ,
    -3.899445819932787e-01,    -1.631991173051828e-01 ,
    -3.889260208856150e-01,     1.002863473536458e+00 ,
    -3.672374225893464e-01,     2.221821776944201e-01 ,
    -3.508504254209158e-01,     5.957889697901096e-01 ,
    -3.213953123774163e-01,    -3.512498585762419e-01 ,
    -2.805213224129314e-01,     2.667423111463128e-02 ,
    -2.761648786285214e-01,    -7.150539086468221e-01 ,
    -2.551406626362323e-01,     3.998534464886444e-01 ,
    -2.494287246067368e-01,    -1.001847600671415e+00 ,
    -2.250077667486312e-01,     7.737609847528025e-01 ,
    -1.921969269222660e-01,    -1.658073724266970e-01 ,
    -1.903379448114237e-01,    -5.375363972704199e-01 ,
    -1.587084461247093e-01,     2.111435963113580e-01 ,
    -1.273356959845733e-01,     5.640562008650558e-01 ,
    -1.078443260647736e-01,    -3.548802371778556e-01 ,
    -7.377816776769955e-02,     2.213039253609772e-02 ,
    -4.654812292827808e-02,    -7.625411824164996e-01 ,
    -3.274041195915585e-02,     3.850544732948923e-01 ,
     7.415424759958738e-03,    -5.552435620481777e-01 ,
     7.682435312354997e-03,    -1.673115638763533e-01 ,
     1.359567110235472e-02,     7.300920540751936e-01 ,
     3.116321947973877e-02,     1.001932017968959e+00 ,
     5.444307805136447e-02,     1.988015103266816e-01 ,
     9.242748289437659e-02,    -3.562405130795229e-01 ,
     1.082588894967674e-01,     5.541224234182670e-01 ,
     1.357445078500805e-01,     1.029243181818670e-02 ,
     1.697111601093526e-01,    -1.001528018598704e+00 ,
     1.847107643500637e-01,     3.686766606118979e-01 ,
     1.965989548482610e-01,    -7.224191701609871e-01 ,
     2.033030530198447e-01,    -5.309632767048605e-01 ,
     2.099633862905836e-01,    -1.765919489160048e-01 ,
     2.573186363992939e-01,     7.601535286062874e-01 ,
     2.681675997008618e-01,     1.727277321402180e-01 ,
     3.031488977338145e-01,    -3.589722788393032e-01 ,
     3.195757781651303e-01,     5.512158667911040e-01 ,
     3.517100700748270e-01,    -1.468397341694293e-02 ,
     3.866386939311920e-01,     3.338701111544212e-01 ,
     3.947605857257251e-01,    -5.516063305296033e-01 ,
     4.194008437796602e-01,    -1.971004466239816e-01 ,
     4.284719840347163e-01,    -7.593392764281580e-01 ,
     4.463830555824670e-01,     1.000096990692893e+00 ,
     4.865237129438513e-01,     1.253285124204812e-01 ,
     5.076515167964308e-01,     5.054848668073821e-01 ,
     5.126747062688567e-01,    -3.679101471246932e-01 ,
     5.127894871146351e-01,     7.135581559853690e-01 ,
     5.715523984380941e-01,    -9.967712133649619e-01 ,
     5.737963184990295e-01,    -9.669118790251183e-02 ,
     5.870002032045341e-01,    -5.681964808138850e-01 ,
     5.889433174409875e-01,     3.129157411249104e-01 ,
     6.919361765873957e-01,    -3.217298746683477e-01 ,
     6.962161115076618e-01,     1.018570156873180e-01 ,
     6.966344876321213e-01,     5.415755536429060e-01 ,
     6.988949780654984e-01,    -7.445820595221496e-01 ,
     7.316871062897004e-01,     7.330978209871578e-01 ,
     7.875892439616202e-01,    -5.363003737046200e-01 ,
     7.877594979450218e-01,     9.571011771205623e-01 ,
     7.919729853614763e-01,    -1.150326213359235e-01 ,
     7.979868908847152e-01,     3.246793059211840e-01 ,
     8.628322308675166e-01,    -9.166691870072856e-01 ,
     9.266580056149254e-01,     8.484084706239443e-01 ,
     9.717072934799963e-01,    -7.423760037989784e-01 ,
     9.973204490577126e-01,     5.463994765975289e-01 ,
     1.002382072710754e+00,    -3.344310858882624e-01 ,
     1.002525609471183e+00,     1.048850020565920e-01 ,
};


/* Triangles */ 
#define NR_DOMAINS 155
#define SIZE_ia 465

static const matlib_index ia[SIZE_ia] =
{
    41,    51,    53, 
    70,    64,    51,
    31,    21,    33,
    85,    77,    83,
    14,    11,     5,
    15,    11,    21,
    14,    33,    21,
    36,    41,    43,
    33,    51,    41,
    53,    64,    62,
    70,    77,    64,
    79,    85,    87,
    72,    64,    77,
    70,    83,    77,
    43,    53,    54,
    51,    64,    53,
     3,    11,     9,
    14,    21,    11,
    24,    31,    36,
    33,    41,    31,
     5,    11,     3,
     9,    15,    19,
    21,    24,    15,
    62,    72,    68,
    77,    79,    72,
    39,    43,    48,
    41,    53,    43,
    62,    68,    58,
    64,    72,    62,
    74,    72,    79,
    77,    85,    79,
    48,    54,    58,
    53,    62,    54,
    29,    36,    39,
    31,    41,    36,
    19,    24,    29,
    21,    31,    24,
     0,     9,     8,
    11,    15,     9,
    44,    39,    48,
    43,    54,    48,
    55,    48,    58,
    54,    62,    58,
    58,    68,    63,
    72,    74,    68,
    35,    29,    39,
    36,    43,    39,
    25,    19,    29,
    24,    36,    29,
     8,    19,    17,
    15,    24,    19,
    74,    81,    71,
    79,    87,    74,
     3,     9,     0,
    55,    63,    60,
    68,    71,    63,
     1,     8,    12,
     9,    19,     8,
    44,    55,    50,
    58,    63,    55,
    35,    44,    40,
    48,    55,    44,
    25,    35,    30,
    39,    44,    35,
    25,    30,    22,
    29,    35,    25,
    12,    17,    22,
    19,    25,    17,
    75,    71,    81,
    74,    87,    81,
    60,    71,    66,
    68,    74,    71,
    50,    60,    57,
    63,    71,    60,
    40,    50,    47,
    55,    60,    50,
    30,    40,    37,
    44,    50,    40,
    22,    30,    27,
    35,    40,    30,
     6,     1,    12,
     8,    17,    12,
    16,    12,    22,
    17,    25,    22,
     0,     8,     1,
    75,    82,    73,
    81,    88,    75,
    81,    87,    88,
    66,    73,    61,
    71,    75,    66,
    47,    57,    52,
    60,    66,    57,
    37,    47,    42,
    50,    57,    47,
    27,    37,    32,
    40,    47,    37,
    16,    27,    23,
    30,    37,    27,
    16,    23,    13,
    22,    27,    16,
     2,     1,     6,
    12,    16,     6,
    67,    61,    73,
    66,    75,    73,
    76,    73,    82,
    75,    88,    82,
    59,    52,    61,
    57,    66,    61,
    42,    52,    49,
    57,    61,    52,
    32,    42,    38,
    47,    52,    42,
    23,    32,    28,
    37,    42,    32,
    13,    23,    18,
    27,    32,    23,
    10,     2,    13,
     6,    16,    13,
    59,    67,    69,
    73,    76,    67,
     6,    13,     2,
    76,    78,    69,
    82,    86,    76,
    78,    86,    84,
    82,    88,    86,
    49,    59,    56,
    61,    67,    59,
    38,    49,    45,
    52,    59,    49,
    28,    38,    34,
    42,    49,    38,
    10,    18,    20,
    23,    28,    18,
    28,    34,    20,
    32,    38,    28,
     4,    10,     7,
    13,    18,    10,
    56,    69,    65,
    67,    76,    69,
    34,    45,    46,
    49,    56,    45,
    80,    69,    78,
    76,    86,    78,
    56,    65,    46,
    59,    69,    56,
     7,    20,    26,
    18,    28,    20,
    34,    46,    26,
    38,    45,    34,
     2,    10,     4,
    78,    84,    80,
    10,    20,     7,
    69,    80,    65,
    45,    56,    46,
    20,    34,    26,
};

#define NR_QNODES 6

static const matlib_real xi4[2*NR_QNODES] = 
{
    -1.081030181680703e-01,    -7.837939636638594e-01, 
    -1.081030181680702e-01,    -1.081030181680701e-01,
    -7.837939636638596e-01,    -1.081030181680699e-01,
    -8.168475729804587e-01,    -8.168475729804582e-01,
     6.336951459609168e-01,    -8.168475729804585e-01,
    -8.168475729804581e-01,     6.336951459609170e-01,
};

static const matlib_real qwa[NR_QNODES] = 
{
     4.467631793560230e-01, 
     4.467631793560230e-01,
     4.467631793560230e-01,
     2.199034873106438e-01,
     2.199034873106438e-01,
     2.199034873106438e-01,
};

#define SIZE_ca 310
static const matlib_real ca[SIZE_ca] = 
{
     1.065873306764452e-01,    -8.288294570587302e-01, 
     3.899118475273877e-01,    -9.192128361306079e-01,
    -3.376681430893215e-01,    -8.265657980345860e-01,
     8.444781674710038e-01,    -8.012090834428044e-01,
    -7.558664285578827e-01,    -8.705348490231889e-01,
    -6.050024834776122e-01,    -6.879826068209921e-01,
    -4.616320097953464e-01,    -9.181222439333802e-01,
    -7.649021432658101e-02,    -6.184403805783657e-01,
    -4.208856247522075e-02,    -9.219722672288727e-01,
     3.399438415362341e-01,    -6.777882590395828e-01,
     5.663064535127696e-01,    -8.335641831050898e-01,
     9.205595367174567e-01,    -5.377024877972870e-01,
     5.714557217682495e-01,    -6.907059389213975e-01,
     7.110932024570364e-01,    -8.860074866314657e-01,
     1.357724775426881e-01,    -6.028753363046752e-01,
     2.649273663307767e-01,    -8.277621550626164e-01,
    -8.243180041877490e-01,    -6.215939323019805e-01,
    -6.228135780947648e-01,    -8.259305605384766e-01,
    -2.920492908010952e-01,    -5.938945994368022e-01,
    -1.907139087211787e-01,    -8.264808972449122e-01,
    -8.681417880958041e-01,    -7.671183456280218e-01,
    -6.328881399446932e-01,    -4.611765976638445e-01,
    -4.972263566303948e-01,    -6.225895874564791e-01,
     4.981451650663719e-01,    -4.959043194893939e-01,
     6.911614750772176e-01,    -6.163596380135515e-01,
    -2.667139470146088e-03,    -4.221214374351854e-01,
     5.248875222664722e-02,    -6.800679715418881e-01,
     4.035280632427988e-01,    -4.261629188311999e-01,
     4.700775909883251e-01,    -6.263806959238821e-01,
     6.888418745845167e-01,    -4.754089097289509e-01,
     8.193971718357050e-01,    -6.744194790085826e-01,
     1.996264778826786e-01,    -4.153920228745622e-01,
     2.648875311979436e-01,    -6.016629257984837e-01,
    -2.065258610845379e-01,    -4.145554976748391e-01,
    -1.710169821227411e-01,    -6.717104961112472e-01,
    -4.266611868606011e-01,    -4.093194161695516e-01,
    -3.910735845415226e-01,    -6.689810952751692e-01,
    -8.389295182627811e-01,    -3.193589752360291e-01,
    -6.942315500326751e-01,    -5.870624819796632e-01,
    -2.578135952680671e-03,    -2.928107713779106e-01,
     1.010486535580600e-01,    -4.808157839441871e-01,
     2.018465889729249e-01,    -2.972682469449436e-01,
     3.004041788264614e-01,    -4.805139620245891e-01,
     4.117414825941104e-01,    -3.079942908626593e-01,
     5.972036953535955e-01,    -4.192788342023086e-01,
    -2.071455217881520e-01,    -2.906458227269315e-01,
    -9.692228203874619e-02,    -4.825533988321511e-01,
    -4.200943645372472e-01,    -2.873546244735576e-01,
    -3.071261020507267e-01,    -4.726265827466089e-01,
    -6.259388937576283e-01,    -2.247625493340493e-01,
    -5.177371476998417e-01,    -4.841959250410550e-01,
     6.859018268159671e-01,    -1.778178946355944e-01,
     8.273024977532565e-01,    -3.974871114204100e-01,
    -9.141925685555087e-01,    -4.829645729470120e-01,
     3.270247667150236e-01,    -1.294587896523098e-01,
     5.019572895158487e-01,    -2.205672605503955e-01,
    -8.229651223428456e-01,    -2.141966622818514e-02,
    -6.877115438070396e-01,    -3.321024502751789e-01,
     1.177967764843397e-01,    -1.112036936580571e-01,
     3.108377092680195e-01,    -2.442215581264299e-01,
    -8.609755312587020e-02,    -1.036628479223175e-01,
     1.033577681657717e-01,    -2.333813419572937e-01,
    -2.875542771094921e-01,    -1.007774195390828e-01,
    -9.745293922489488e-02,    -2.293330578269686e-01,
    -3.852291265278567e-01,    -2.912443088730883e-02,
    -3.011789404309870e-01,    -2.267521161027072e-01,
    -5.876456779175776e-01,    -1.456385806880157e-02,
    -5.028892855946622e-01,    -2.162766074275253e-01,
     6.873284717893892e-01,    -3.662226451703911e-02,
     8.287637448865420e-01,    -2.570645272975112e-01,
     4.706767005059025e-01,     4.651117033675485e-03,
     5.928024004517606e-01,    -2.621104032318509e-01,
     2.518740592085898e-01,     5.611206351382059e-02,
     4.483024107845056e-01,    -1.028252026478121e-01,
     3.880313937791514e-02,     7.707477822698867e-02,
     2.324726547384971e-01,    -6.032783017158701e-02,
    -1.710026454351134e-01,     8.664940665402900e-02,
     2.321625846491198e-02,    -4.496291317402296e-02,
    -3.776600733932127e-01,     9.933600077922548e-02,
    -1.821654723676323e-01,    -3.900091625865600e-02,
    -8.363371262165716e-01,     1.337675531767378e-01,
    -6.789363216858501e-01,    -9.383337025326155e-02,
    -5.920523984046351e-01,     1.127881325910962e-01,
    -4.816487109067666e-01,    -8.402111040490091e-02,
    -9.211856688703653e-01,    -1.396053702698231e-01,
     6.943821066111214e-01,     2.464840209111374e-01,
     8.302382354467737e-01,     3.056979880266217e-02,
     9.322935558478044e-01,    -1.148595683891980e-01,
     4.873685747720103e-01,     2.573714548999376e-01,
     5.855120476501808e-01,     4.349811340176246e-02,
     1.691071473674300e-01,     2.467353010262658e-01,
     3.688004609065134e-01,     9.445742371458543e-02,
    -4.566859334416690e-02,     2.649998599776440e-01,
     1.527850618674356e-01,     1.272738914283621e-01,
    -2.603621771167626e-01,     2.777264068314741e-01,
    -5.934784528034812e-02,     1.440251663913791e-01,
    -4.731104824625765e-01,     2.923171477501467e-01,
    -2.688223970423290e-01,     1.533333350401365e-01,
    -5.851814704173658e-01,     3.736975810205842e-01,
    -4.784863782591800e-01,     1.717913119215312e-01,
    -9.331805078587704e-01,     2.933085574895355e-01,
    -6.967150460466335e-01,     1.887058398115592e-01,
     4.944111760562035e-01,     3.840902396955713e-01,
     5.905610472975001e-01,     1.800337564109032e-01,
     6.945215653192746e-01,     3.930568668963335e-01,
     8.322428706211866e-01,     1.771404412216980e-01,
     2.969750788154620e-01,     4.179208795191410e-01,
     3.804433355253017e-01,     2.106421185717068e-01,
     8.674308062922509e-02,     4.359511857750191e-01,
     2.798390193273725e-01,     2.917581679688457e-01,
    -1.384055901933205e-01,     4.496547068828642e-01,
     6.880447681409078e-02,     3.175108814111573e-01,
    -3.583616252815658e-01,     4.687905057644086e-01,
    -1.488631735733658e-01,     3.320171720316316e-01,
    -5.788764550236943e-01,     4.876669760907335e-01,
    -3.638239576710428e-01,     3.442549083991787e-01,
    -8.162570320668546e-01,     5.653448550527060e-01,
    -6.952200138559678e-01,     3.290894524124317e-01,
     4.466722606920654e-01,     5.900862965279517e-01,
     5.977431072898466e-01,     4.533253871917329e-01,
    -8.337084775031257e-01,     4.223412817790851e-01,
     6.470370270121523e-01,     6.627438435384776e-01,
     8.306472758581830e-01,     4.708847787205396e-01,
     8.852218536541128e-01,     7.093019227362104e-01,
     9.326109831378702e-01,     3.253212615251016e-01,
     2.283844346870639e-01,     6.218306062718861e-01,
     4.046219962975844e-01,     4.635236149176358e-01,
    -1.827045128483729e-03,     6.160902261195055e-01,
     2.041818106706538e-01,     4.913383169404229e-01,
    -2.343979627180401e-01,     6.445353851359893e-01,
    -1.727240614898725e-02,     5.010776991927384e-01,
    -6.111589715997998e-01,     6.881065660324288e-01,
    -4.613431346794283e-01,     5.308221401855258e-01,
    -3.661316114481951e-01,     7.140019559591085e-01,
    -2.444422613472405e-01,     5.198995390479366e-01,
    -8.141900202153561e-01,     8.391434497828104e-01,
    -6.714635530260251e-01,     5.860624207562019e-01,
     4.054970596987987e-01,     8.246028917615164e-01,
     5.723584971810625e-01,     5.868728588118857e-01,
    -6.008295872217925e-02,     8.352616855989851e-01,
     1.263910656661387e-01,     6.814560020332493e-01,
     6.774120304497858e-01,     8.012523846976963e-01,
     8.085473476598448e-01,     6.070242837425308e-01,
     2.449549704871666e-01,     9.207275124227131e-01,
     3.632279672263531e-01,     6.749758504609201e-01,
    -5.601148305169786e-01,     9.128653652656172e-01,
    -4.791574194752579e-01,     6.513977442921730e-01,
    -1.942568560515025e-01,     9.261854920860731e-01,
    -1.129159305436166e-01,     6.893030798976839e-01,
    -8.907179533693647e-01,     6.926497643546533e-01,
     8.153682032832158e-01,     8.462024895772214e-01,
    -6.794245174866207e-01,     8.138826857570901e-01,
     5.823106802140413e-01,     8.902521079329414e-01,
     1.006925089937958e-01,     8.307258668834802e-01,
    -3.788234766030948e-01,     8.496934572078914e-01,
};



static const fem2d_cc nodes = { .len = SIZE_p/2, 
                                .elem_p = (matlib_real*) p };

static const fem2d_cc cen = { .len = NR_DOMAINS, 
                              .elem_p = (matlib_real*)ca };

static const fem2d_cc xi = { .len = NR_QNODES, 
                             .elem_p = (matlib_real*)xi4 };

/*============================================================================*/

int init_suite(void)
{
      return 0;
}
int clean_suite(void)
{
      return 0;
}              
/*============================================================================*/
static fem2d_err cyclic_permutation
(
    const matlib_index *const ia, 
    const matlib_index len,
    const matlib_index p_step,
         matlib_index* cia
)
{
    int i;
    fem2d_err error = FEM2D_SUCCESS;

    if (p_step == 0)
    {
        for (i = 0; i < 3*len; i+=3)
        {
            cia[0 + i] = ia[0 + i];
            cia[1 + i] = ia[1 + i];
            cia[2 + i] = ia[2 + i];
        }
    }
    if (p_step == 1)
    {
        for (i = 0; i < 3*len; i+=3)
        {
            cia[0 + i] = ia[2 + i];
            cia[1 + i] = ia[0 + i];
            cia[2 + i] = ia[1 + i];
        }
    }
    else if (p_step == 2)
    {
        for (i = 0; i < 3*len; i+=3)
        {
            cia[0 + i] = ia[1 + i];
            cia[1 + i] = ia[2 + i];
            cia[2 + i] = ia[0 + i];
        }
    }
    else
    {
        error = FEM2D_FAILURE;
    }
    return error;
}
/*============================================================================*/

void test_create_ea(void)
{
    fem2d_ea ea;
    fem2d_create_ea(nodes, ia, NR_DOMAINS, &ea);


    fem2d_te* dptr = ea.elem_p; /* domain pointer */
    matlib_real* vptr; /* vertex pointer */ 

    matlib_index const *pia = ia;
    matlib_real e_relative, x1, xx1, x2, xx2;

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); 
          dptr++, pia +=3)
    {
        /* Get the first vertex */ 
        vptr = *(dptr->vert_p);
        
        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *pia + 0];
        xx2 = p[2 * *pia + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);

        /* Get the second vertex */
        vptr = *(dptr->vert_p+1);

        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *(pia+1) + 0];
        xx2 = p[2 * *(pia+1) + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);
        
        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);

        /* Get the third vertex */
        vptr = *(dptr->vert_p+2);

        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *(pia+2) + 0];
        xx2 = p[2 * *(pia+2) + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}

void test_create_ea1(void)
{
    matlib_index cia[SIZE_ia];
    cyclic_permutation(ia, NR_DOMAINS, 1, cia);

    fem2d_ea ea;
    fem2d_create_ea(nodes, cia, NR_DOMAINS, &ea);
   
    fem2d_te* dptr = ea.elem_p; /* domain pointer */
    matlib_real* vptr; /* vertex pointer */ 

    matlib_index *pia = cia;
    matlib_real e_relative, x1, xx1, x2, xx2;

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); 
          dptr++, pia +=3)
    {
        /* Get the first vertex */ 
        vptr = *(dptr->vert_p+0);
        
        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *pia + 0];
        xx2 = p[2 * *pia + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);

        /* Get the second vertex */
        vptr = *(dptr->vert_p+1);

        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *(pia+1) + 0];
        xx2 = p[2 * *(pia+1) + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);
        
        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);

        /* Get the third vertex */
        vptr = *(dptr->vert_p+2);

        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *(pia+2) + 0];
        xx2 = p[2 * *(pia+2) + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}

void test_create_ea2(void)
{
    matlib_index cia[SIZE_ia];
    cyclic_permutation(ia, NR_DOMAINS, 2, cia);

    fem2d_ea ea;
    fem2d_create_ea(nodes, cia, NR_DOMAINS, &ea);
   
    fem2d_te* dptr = ea.elem_p; /* domain pointer */
    matlib_real* vptr; /* vertex pointer */ 

    matlib_index *pia = cia;
    matlib_real e_relative, x1, xx1, x2, xx2;

    for ( dptr = ea.elem_p; 
          dptr < (ea.len + ea.elem_p); 
          dptr++, pia +=3)
    {
        /* Get the first vertex */ 
        vptr = *(dptr->vert_p+0);
        
        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *pia + 0];
        xx2 = p[2 * *pia + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);

        /* Get the second vertex */
        vptr = *(dptr->vert_p+1);

        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *(pia+1) + 0];
        xx2 = p[2 * *(pia+1) + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);
        
        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);

        /* Get the third vertex */
        vptr = *(dptr->vert_p+2);

        x1 = vptr[0]; 
        x2 = vptr[1]; 
        
        xx1 = p[2 * *(pia+2) + 0];
        xx2 = p[2 * *(pia+2) + 1];
        
        e_relative = fabs(x1-xx1)/(TOL+xx1);
        CU_ASSERT_TRUE(e_relative<TOL);

        e_relative = fabs(x2-xx2)/(TOL+xx2);
        CU_ASSERT_TRUE(e_relative<TOL);
    }
}
/*============================================================================*/
void test_create_ia(void)
{

    fem2d_ea ea;
    fem2d_create_ea(nodes, ia, NR_DOMAINS, &ea);

    matlib_index ia_tmp[SIZE_ia];
    fem2d_create_ia(ea, ia_tmp);
    
    BEGIN_DEBUG
        for (matlib_index i=0; i< SIZE_ia; i++)
        {
            debug_print("ia[%d] -> actual: %d, computed: %d", i, ia[i], ia_tmp[i]);
        }
    END_DEBUG

    CU_ASSERT_TRUE(true);
}


/*============================================================================*/
void test_centroid(void)
{
    

    fem2d_ea ea0;
    fem2d_create_ea(nodes, ia, NR_DOMAINS, &ea0);
    
    fem2d_cc cen0;
    fem2d_create_cc(NR_DOMAINS, &cen0);

    fem2d_centroid(ea0, cen0);
    
    matlib_real err;
    int i;
    for (i = 0; i< 2*NR_DOMAINS; i+=2)
    {
        err =   fabs(cen0.elem_p[i]-cen.elem_p[i]) 
              + fabs(cen0.elem_p[i+1]-cen.elem_p[i+1]);
              
        CU_ASSERT_TRUE(err<TOL);
    }
}

void test_centroid1(void)
{

    matlib_index cia1[SIZE_ia], cia2[SIZE_ia];
    cyclic_permutation(ia, NR_DOMAINS, 1, cia1);
    cyclic_permutation(ia, NR_DOMAINS, 2, cia2);

    fem2d_ea ea0, ea1, ea2;
    fem2d_create_ea(nodes,   ia, NR_DOMAINS, &ea0);
    fem2d_create_ea(nodes, cia1, NR_DOMAINS, &ea1);
    fem2d_create_ea(nodes, cia2, NR_DOMAINS, &ea2);
    
    fem2d_cc cen0, cen1, cen2;
    fem2d_create_cc(NR_DOMAINS, &cen0);
    fem2d_create_cc(NR_DOMAINS, &cen1);
    fem2d_create_cc(NR_DOMAINS, &cen2);

    fem2d_centroid(ea0, cen0);
    fem2d_centroid(ea1, cen1);
    fem2d_centroid(ea2, cen2);
    
    matlib_real err;
    int i;
    for (i = 0; i< 2*NR_DOMAINS; i+=2)
    {
        err = fabs(cen0.elem_p[i]-cen1.elem_p[i])     + fabs(cen1.elem_p[i]-cen2.elem_p[i]) + 
              fabs(cen0.elem_p[i+1]-cen1.elem_p[i+1]) + fabs(cen1.elem_p[i+1]-cen2.elem_p[i+1]); 
              
        CU_ASSERT_TRUE(err<TOL);
    }
}


/*============================================================================*/
void test_refbasis(void)
{

    matlib_xm vphi;
    

    fem2d_refbasis(xi, &vphi);
    
    int i;
    
    matlib_real err;

    for(i = 0; i < xi.len; i++)
    {
        err =   vphi.elem_p[i] 
              + vphi.elem_p[xi.len + i] 
              + vphi.elem_p[2*xi.len + i] 
              - 1.0;
        
        CU_ASSERT_TRUE(fabs(err)<TOL);
    }
}


void test_ref2mesh(void)
{
    matlib_index cia1[SIZE_ia], cia2[SIZE_ia];
    cyclic_permutation(ia, NR_DOMAINS, 1, cia1);
    cyclic_permutation(ia, NR_DOMAINS, 2, cia2);

    fem2d_ea ea0, ea1, ea2;
    fem2d_create_ea(nodes,   ia, NR_DOMAINS, &ea0);
    fem2d_create_ea(nodes, cia1, NR_DOMAINS, &ea1);
    fem2d_create_ea(nodes, cia2, NR_DOMAINS, &ea2);

    matlib_xm vphi;
    fem2d_refbasis(xi, &vphi);

    fem2d_cc x0, x1, x2;
    fem2d_ref2mesh(ea0, vphi, &x0);
    fem2d_ref2mesh(ea1, vphi, &x1);
    fem2d_ref2mesh(ea2, vphi, &x2);
    
    matlib_index i, j; 
    matlib_real err1, err2, cent1, cent2;

    for (i = 0; i < ea0.len; i++ )
    {
        cent1 = ca[2 * i + 0];
        cent2 = ca[2 * i + 1];

        for (j = 0; j < (2 * xi.len); j+=2 )
        {
            err1 =  fabs(  x0.elem_p[2 * i * xi.len + j + 0] 
                         + x1.elem_p[2 * i * xi.len + j + 0] 
                         + x2.elem_p[2 * i * xi.len + j + 0]
                         - 3.0 *  cent1); 
            err2 =  fabs(  x0.elem_p[2 * i * xi.len + j + 1] 
                         + x1.elem_p[2 * i * xi.len + j + 1] 
                         + x2.elem_p[2 * i * xi.len + j + 1]
                         - 3.0 *  cent2); 

            debug_body("domain: %d, error: %0.16f", i,  err1 + err2);

            CU_ASSERT_TRUE((err1 + err2) < TOL);
        }   
    }
}

/*============================================================================*/
fem2d_err poly_func
(
    fem2d_cc nodes_tmp, 
    matlib_zv u_nodes
)
{
    debug_enter( "nodes: %d, u_nodes: %d", 
                 nodes_tmp.len, u_nodes.len);
    err_check(( nodes_tmp.len != u_nodes.len), clean_up, 
                "Dimension mismatch (nodes: %d, u_nodes: %d)!", 
                nodes_tmp.len, u_nodes.len);
    
    matlib_real* ptr;
    matlib_complex* uptr = u_nodes.elem_p;
    for ( ptr = nodes_tmp.elem_p; 
          ptr < (nodes_tmp.elem_p +2 * (nodes_tmp.len)); ptr+=2, uptr++)
    {
        *uptr = (ptr[0] + ptr[1]);
    }
    debug_exit("exit status: %d", FEM2D_SUCCESS);
    return FEM2D_SUCCESS;
clean_up:
    debug_exit("exit status: %d", FEM2D_FAILURE);
    return FEM2D_FAILURE;
}

void test_interp(void)
{

    fem2d_ea ea;
    fem2d_create_ea(nodes, ia, NR_DOMAINS, &ea);

    matlib_zv u_nodes;
    matlib_create_zv( nodes.len, &u_nodes, MATLIB_COL_VECT);
    poly_func(nodes, u_nodes);
    DEBUG_PRINT_ZV(u_nodes, "%s:", "u_nodes");

    matlib_xm vphi;
    fem2d_refbasis(xi, &vphi);
    DEBUG_PRINT_XM( vphi, "%s:", "vphi");

    fem2d_cc x_interp;
    fem2d_ref2mesh(ea, vphi, &x_interp);

    matlib_zv u_interp;
    matlib_zv u_exact;

    matlib_create_zv( NR_DOMAINS * NR_QNODES, &u_interp, MATLIB_COL_VECT);
    matlib_create_zv( NR_DOMAINS * NR_QNODES, &u_exact , MATLIB_COL_VECT);
    BEGIN_DEBUG
        for (matlib_index i = 0; i < u_interp.len; i++)
        {
            debug_print( "u[%d]-> diff: %0.16f%+0.16fi",
                         i, 
                         creal(u_interp.elem_p[i])-creal(u_exact.elem_p[i]),
                         cimag(u_interp.elem_p[i])-cimag(u_exact.elem_p[i]));
        
        }
    END_DEBUG

    fem2d_interp( ea, u_nodes, vphi, u_interp);
    poly_func(x_interp, u_exact);

 
    matlib_real norm_actual = matlib_znrm2(u_exact);
    debug_body("norm: %0.16g", norm_actual);
    matlib_zaxpy(-1.0, u_exact, u_interp);
    matlib_real e_relative = matlib_znrm2(u_interp)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    CU_ASSERT_TRUE(e_relative<TOL);

}

fem2d_err test_interp1(void)
{

    fem2d_err error = FEM2D_SUCCESS;

    matlib_index mcnt_total = 6;
    matlib_index mcnt = 0;

    void** ptr = calloc(mcnt_total, sizeof(void*));

    /* ea: element array */ 
    ptr[mcnt] = calloc(1, sizeof(fem2d_ea));
    fem2d_ea ea = *((fem2d_ea*)ptr[0]);
    error = fem2d_create_ea(nodes, ia, NR_DOMAINS, &ea);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Creation of element array failed!");
    mcnt++;

    /* u_nodes: field values at nodes */ 
    ptr[mcnt] = calloc(1, sizeof(matlib_zv));
    matlib_zv u_nodes = *((matlib_zv*)ptr[mcnt]);
    matlib_create_zv( nodes.len, &u_nodes, MATLIB_COL_VECT);
    mcnt++;
    error = poly_func(nodes, u_nodes);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Calculation of field at nodes failed!");
    DEBUG_PRINT_ZV(u_nodes, "%s:", "u_nodes");

    /* vphi: matrix containing values of reference basis functions at 
     * quadrature points
     * */ 
    ptr[mcnt] = calloc(1, sizeof(matlib_zv));
    matlib_xm vphi = *((matlib_xm*)ptr[mcnt]);
    error = fem2d_refbasis(xi, &vphi);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Creation of basis function-value matrix failed!");
    mcnt++;
    DEBUG_PRINT_XM( vphi, "%s:", "vphi");

    /* x_interp: mesh points for interpolation 
     * */ 
    ptr[mcnt] = calloc(1, sizeof(fem2d_cc));
    fem2d_cc x_interp = *((fem2d_cc*)ptr[mcnt]);
    error = fem2d_ref2mesh(ea, vphi, &x_interp);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Creation of interpolation mesh-point array failed!");
    mcnt++;

    ptr[mcnt] = calloc(1, sizeof(matlib_zv));
    matlib_zv u_interp = *((matlib_zv*)ptr[mcnt]);
    mcnt++;

    ptr[mcnt] = calloc(1, sizeof(matlib_zv));
    matlib_zv u_exact = *((matlib_zv*)ptr[mcnt]);
    mcnt++;

    matlib_create_zv( NR_DOMAINS * NR_QNODES, &u_interp, MATLIB_COL_VECT);
    matlib_create_zv( NR_DOMAINS * NR_QNODES, &u_exact , MATLIB_COL_VECT);
    BEGIN_DEBUG
        for (matlib_index i = 0; i < u_interp.len; i++)
        {
            debug_print( "u[%d]-> diff: %0.16f%+0.16fi",
                         i, 
                         creal(u_interp.elem_p[i])-creal(u_exact.elem_p[i]),
                         cimag(u_interp.elem_p[i])-cimag(u_exact.elem_p[i]));
        
        }
    END_DEBUG

    error = fem2d_interp( ea, u_nodes, vphi, u_interp);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Interpolation failed!");

    error = poly_func(x_interp, u_exact);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Calculation of exact field values failed!");
 
    matlib_real norm_actual = matlib_znrm2(u_exact);
    matlib_zaxpy(-1.0, u_exact, u_interp);
    matlib_real e_relative = matlib_znrm2(u_interp)/norm_actual;

    debug_exit("Relative error: %0.16g", e_relative);
    err_check( (e_relative > TOL), clean_up, "%s", 
               "Tolerance test failed!");
    
    
    matlib_free(u_exact.elem_p);
    matlib_free(u_interp.elem_p);
    matlib_free(ptr[5]);
    matlib_free(ptr[4]);
    matlib_free(x_interp.elem_p);
    matlib_free(ptr[3]);
    matlib_free(vphi.elem_p);
    matlib_free(ptr[2]);
    matlib_free(u_nodes.elem_p);
    matlib_free(ptr[1]);
    fem2d_free_ea(ea);
    matlib_free(ptr[0]);
    matlib_free(ptr);
    return FEM2D_SUCCESS;

clean_up:
    if(mcnt==5)
    {
        matlib_free(u_exact.elem_p);
        matlib_free(u_interp.elem_p);
        matlib_free(ptr[5]);
        matlib_free(ptr[4]);
        mcnt -= 2;
    }
    if(mcnt==3)
    {
        matlib_free(x_interp.elem_p);
        matlib_free(ptr[3]);
        mcnt--;
    }
    if(mcnt==2)
    {
        matlib_free(vphi.elem_p);
        matlib_free(ptr[2]);
        mcnt--;
    }
    if(mcnt==1)
    {
        matlib_free(u_nodes.elem_p);
        matlib_free(ptr[1]);
        mcnt--;
    }
    if(mcnt==0)
    {
        fem2d_free_ea(ea);
        matlib_free(ptr[0]);
    }
    matlib_free(ptr);
    return FEM2D_FAILURE;
}

void test_normL2(void)
{

    fem2d_err error = FEM2D_SUCCESS;

    matlib_index mcnt_total = 6;
    matlib_index mcnt = 0;

    void** ptr = calloc(mcnt_total, sizeof(void*));

    /* ea: element array */ 
    errno = 0;
    ptr[mcnt] = calloc(1, sizeof(fem2d_ea));
    fem2d_ea ea = *((fem2d_ea*)ptr[0]);
    error = fem2d_create_ea(nodes, ia, NR_DOMAINS, &ea);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Creation of element array failed!");
    mcnt++; /* 1 */ 

    /* u_nodes: field values at nodes */ 
    ptr[mcnt] = calloc(1, sizeof(matlib_zv));
    matlib_zv u_nodes = *((matlib_zv*)ptr[mcnt]);
    matlib_create_zv( nodes.len, &u_nodes, MATLIB_COL_VECT);
    mcnt++; /* 2 */ 
    error = poly_func(nodes, u_nodes);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Calculation of field at nodes failed!");
    DEBUG_PRINT_ZV(u_nodes, "%s:", "u_nodes");

    /* vphi: matrix containing values of reference basis functions at 
     * quadrature points
     * */ 
    ptr[mcnt] = calloc(1, sizeof(matlib_zv));
    matlib_xm vphi = *((matlib_xm*)ptr[mcnt]);
    error = fem2d_refbasis(xi, &vphi);
    err_check( (error == FEM2D_FAILURE), clean_up, "%s", 
               "Creation of basis function-value matrix failed!");
    mcnt++; /* 3 */ 
    DEBUG_PRINT_XM( vphi, "%s:", "vphi");

    /* quadW: quadrature weights */ 
    ptr[mcnt] = calloc(1, sizeof(matlib_xv));
    matlib_xv quadW = *((matlib_xm*)ptr[mcnt]);

    mcnt++; /* 3 */ 
    DEBUG_PRINT_XM( vphi, "%s:", "vphi");



}


void test_interp11(void)
{
    fem2d_err error = test_interp1();
    CU_ASSERT_TRUE(error == FEM2D_SUCCESS);

}
/*============================================================================*/

int main()
{
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
    {
        return CU_get_error();
    }

    /* Create a test array */
    CU_TestInfo test_array[] = 
    {
        { "Create element array"               , test_create_ea },
        { "Create element array1"              , test_create_ea1},
        { "Create element array2"              , test_create_ea2},
        { "Create index array"              , test_create_ia},
        { "Centroid"                           , test_centroid  },
        { "Centroid1"                          , test_centroid1 },
        { "Ref. basis"                         , test_refbasis  },
        {"Ref to mesh"                         , test_ref2mesh  },
        { "Interpolation on trianglar elements", test_interp    },
        { "Interpolation on trianglar elements1", test_interp11    },
        CU_TEST_INFO_NULL,
    };

    /* Create the test suite */ 
    CU_SuiteInfo suites[] = 
    {
        { "fem2d", init_suite, clean_suite, NULL, NULL, (CU_TestInfo*) test_array },
        CU_SUITE_INFO_NULL
    }; 

    /* Register test suites */ 
    CU_ErrorCode CU_error = CU_register_suites(suites); 
    if (CU_error != CUE_SUCCESS) 
    {
        debug_body("%s", CU_get_error_msg());
        CU_cleanup_registry();
        return CU_get_error();
    }


   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}
