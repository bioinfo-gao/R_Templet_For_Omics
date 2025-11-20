# CD45(In115)Dd 是 CyTOF 实验中标记抗体的标准命名格式
# 靶点 (CD45)： 目标生物分子。
# 金属标签 (In115)： 用于检测的重金属同位素 铟（Indium） 元素。
# DNA 偶联物 (Dd)： 抗体偶联的技术细节。
# 在 CyTOF 中，这个特定的金属同位素被共价偶联到 CD45 标记物。 

# # \[\d+\] to match [9] or [81]

# "CD3(110:114)Dd"  "CD45(In115)Dd"  
# "BC1(La139)Dd"  "BC2(Pr141)Dd"    "pNFkB(Nd142)Dd"  "pp38(Nd144)Dd"  
# "CD4(Nd145)Dd"    "BC3(Nd146)Dd"    "CD20(Sm147)Dd"   "CD33(Nd148)Dd"  
# "pStat5(Nd150)Dd" "CD123(Eu151)Dd"  "pAkt(Sm152)Dd"   "pStat1(Eu153)Dd"
# "pSHP2(Sm154)Dd"  "pZap70(Gd156)Dd" "pStat3(Gd158)Dd" "BC4(Tb159)Dd"   
# "CD14(Gd160)Dd"   "pSlp76(Dy164)Dd" "BC5(Ho165)Dd"    "pBtk(Er166)Dd"  
# "pPlcg2(Er167)Dd" "pErk(Er168)Dd"   "BC6(Tm169)Dd"    "pLat(Er170)Dd"  
# "IgM(Yb171)Dd"    "pS6(Yb172)Dd"    "HLA-DR(Yb174)Dd" "BC7(Lu175)Dd"   
# "CD7(Yb176)Dd"    "DNA-1(Ir191)Dd"  "DNA-2(Ir193)Dd"  "group_id"       
# "patient_id"      "sample_id"       "population_id"  

# 这个面板总共使用了 35 个用于标记生物标志物的通道，以及 5 个用于样本识别的 ID 通道。A. 表面标志物 (Surface Markers)这些通道用于识别细胞的类型和亚群（主要是淋巴细胞和髓系细胞）：通道标签目标分子生物学作用$\text{CD}3(110:114)\text{Dd}$$\text{Cd}/\text{In}$$\text{CD}3$标记所有 T 淋巴细胞$\text{CD}45(\text{In}115)\text{Dd}$$\text{In}115$$\text{CD}45$标记所有 白细胞$\text{CD}4(\text{Nd}145)\text{Dd}$$\text{Nd}145$$\text{CD}4$标记辅助性 $\text{T}$ 细胞$\text{CD}20(\text{Sm}147)\text{Dd}$$\text{Sm}147$$\text{CD}20$标记 B 淋巴细胞$\text{CD}33(\text{Nd}148)\text{Dd}$$\text{Nd}148$$\text{CD}33$标记 髓系细胞（如单核/粒细胞）$\text{CD}123(\text{Eu}151)\text{Dd}$$\text{Eu}151$$\text{CD}123$标记浆细胞样树突状细胞（$\text{pDC}$）$\text{CD}14(\text{Gd}160)\text{Dd}$$\text{Gd}160$$\text{CD}14$标记 单核细胞$\text{HLA-DR}(\text{Yb}174)\text{Dd}$$\text{Yb}174$$\text{HLA-DR}$标记 $\text{MHC}$ II 类分子，在 $\text{B}$ 细胞、单核细胞和树突状细胞上表达$\text{CD}7(\text{Yb}176)\text{Dd}$$\text{Yb}176$$\text{CD}7$标记 $\text{T}$ 细胞和 $\text{NK}$ 细胞$\text{IgM}(\text{Yb}171)\text{Dd}$$\text{Yb}171$$\text{IgM}$标记 $\text{B}$ 细胞受体B. 细胞内信号通路（磷酸化蛋白）这些通道用于评估细胞的功能状态和信号传导，通常以 $\text{p}$ 开头，表示磷酸化 (Phosphorylation) 状态，这是许多细胞信号通路激活的标志：通道标签目标分子生物学作用$\text{pNFkB}(\text{Nd}142)\text{Dd}$$\text{Nd}142$磷酸化 $\text{NF}-\kappa\text{B}$炎症、免疫反应和细胞存活信号通路$\text{pp}38(\text{Nd}144)\text{Dd}$$\text{Nd}144$磷酸化 $\text{p}38$应激反应和细胞凋亡信号通路$\text{pStat}5(\text{Nd}150)\text{Dd}$$\text{Nd}150$磷酸化 $\text{Stat}5$细胞增殖、分化、生存（如对 IL-7, IL-2 的反应）$\text{pAkt}(\text{Sm}152)\text{Dd}$$\text{Sm}152$磷酸化 $\text{Akt}$细胞生长、增殖和代谢的关键通路 ($\text{PI}3\text{K}/\text{Akt}$)$\text{pStat}1(\text{Eu}153)\text{Dd}$$\text{Eu}153$磷酸化 $\text{Stat}1$干扰素信号和抗病毒反应$\text{pSHP}2(\text{Sm}154)\text{Dd}$$\text{Sm}154$磷酸化 $\text{SHP}2$细胞因子和生长因子受体下游信号$\text{pZap}70(\text{Gd}156)\text{Dd}$$\text{Gd}156$磷酸化 $\text{Zap}70$$\text{T}$ 细胞受体 ($\text{TCR}$) 信号通路的关键$\text{pStat}3(\text{Gd}158)\text{Dd}$$\text{Gd}158$磷酸化 $\text{Stat}3$炎症、增殖和肿瘤信号通路$\text{pSlp}76(\text{Dy}164)\text{Dd}$$\text{Dy}164$磷酸化 $\text{SLP}76$$\text{TCR}$ 和 $\text{BCR}$ 信号传导$\text{pBtk}(\text{Er}166)\text{Dd}$$\text{Er}166$磷酸化 $\text{Btk}$$\text{B}$ 细胞受体 ($\text{BCR}$) 信号通路$\text{pPlcg}2(\text{Er}167)\text{Dd}$$\text{Er}167$磷酸化 $\text{Plc}\gamma 2$$\text{BCR}$ 和其他酪氨酸激酶受体的下游信号$\text{pErk}(\text{Er}168)\text{Dd}$$\text{Er}168$磷酸化 $\text{Erk}$$\text{MAPK}$ 通路，调节增殖和分化$\text{pLat}(\text{Er}170)\text{Dd}$$\text{Er}170$磷酸化 $\text{LAT}$$\text{T}$ 细胞激活的早期信号$\text{pS}6(\text{Yb}172)\text{Dd}$$\text{Yb}172$磷酸化 $\text{S}6$$\text{mTOR}$ 信号通路，调控蛋白质合成和细胞生长C. 细胞和样本识别 ($\text{BC}$ 和 $\text{DNA}$)$\text{BC}1$ 到 $\text{BC}7$ (Barcoding 通道)：$\text{BC}$ 代表 Barcoding（条形码/样本标记）。在 CyTOF 实验中，为了节省时间和试剂，通常会将来自不同条件（如不同药物处理、不同病人）的多个样本混合在一起进行染色和运行。这些 $\text{BC}$ 通道（$\text{La}139, \text{Pr}141$ 等）用于标记样本特有的组合。分析时，数据通过这些 $\text{BC}$ 通道的值被“解复用”，从而区分出原始的每个样本。$\text{DNA}-1(\text{Ir}191)\text{Dd}$ 和 $\text{DNA}-2(\text{Ir}193)\text{Dd}$：$\text{Ir}$ 是铱元素。这两个通道用于标记细胞的 DNA。作用： 主要用于排除细胞碎片和死亡细胞（死亡细胞的 DNA 染色通透性会增强，导致 DNA 信号异常）。D. 数据标识符 (ID)$\text{group\_id}$$\text{patient\_id}$$\text{sample\_id}$$\text{population\_id}$这些是添加到数据中的元数据列，用于在解复用和下游统计分析时追踪每个细胞的来源。

marker_info
channel_name   marker_name marker_class
1             Time          Time         none
2      Cell_length   Cell_length         none
3   CD3(110:114)Dd           CD3         type
4    CD45(In115)Dd          CD45         type
5     BC1(La139)Dd           BC1         none
6     BC2(Pr141)Dd           BC2         none
7   pNFkB(Nd142)Dd         pNFkB        state
8    pp38(Nd144)Dd          pp38        state
9     CD4(Nd145)Dd           CD4         type
10    BC3(Nd146)Dd           BC3         none
11   CD20(Sm147)Dd          CD20         type
12   CD33(Nd148)Dd          CD33         type
13 pStat5(Nd150)Dd        pStat5        state
14  CD123(Eu151)Dd         CD123         type
15   pAkt(Sm152)Dd          pAkt        state
16 pStat1(Eu153)Dd        pStat1        state
17  pSHP2(Sm154)Dd         pSHP2        state
18 pZap70(Gd156)Dd        pZap70        state
19 pStat3(Gd158)Dd        pStat3        state
20    BC4(Tb159)Dd           BC4         none
21   CD14(Gd160)Dd          CD14         type
22 pSlp76(Dy164)Dd        pSlp76        state
23    BC5(Ho165)Dd           BC5         none
24   pBtk(Er166)Dd          pBtk        state
25 pPlcg2(Er167)Dd        pPlcg2        state
26   pErk(Er168)Dd          pErk        state
27    BC6(Tm169)Dd           BC6         none
28   pLat(Er170)Dd          pLat        state
29    IgM(Yb171)Dd           IgM         type
30    pS6(Yb172)Dd           pS6        state
31 HLA-DR(Yb174)Dd        HLA-DR         type
32    BC7(Lu175)Dd           BC7         none
33    CD7(Yb176)Dd           CD7         type
34  DNA-1(Ir191)Dd         DNA-1         none
35  DNA-2(Ir193)Dd         DNA-2         none
36        group_id      group_id         none
37      patient_id    patient_id         none
38       sample_id     sample_id         none
39   population_id population_id         none


