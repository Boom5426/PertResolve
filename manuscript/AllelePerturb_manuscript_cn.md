# AllelePerturb reveals detection limits in single-cell perturbation prediction

李波 (Bo Li)^1

^1 佛罗里达大学, 美国佛罗里达州盖恩斯维尔

---

## 摘要

单细胞扰动预测旨在预测细胞群体对遗传或化学干预的转录组响应,是虚拟细胞建模和精准医学的核心技术之一。现有方法和基准在基因或药物粒度上定义扰动,隐含了同一基因的所有突变产生相同细胞后果这一假设。然而在临床遗传学中,同一基因的不同蛋白编码变异往往导致截然不同的转录调控、药物敏感性和预后。

本文提出 AllelePerturb,首个面向蛋白编码变异水平的单细胞扰动预测基准与评价框架。AllelePerturb 整合了 TP53、KRAS、GATA1 和 JAK1 四个基因共 472 个蛋白编码变异(321,043 个单细胞),覆盖 Perturb-seq、碱基编辑和 scSNV-seq 三种实验技术。我们设计了包含 10 个指标的评价协议 AllelePerturb-Eval,将评测分解为方向恢复、变异排序和差异表达保真度三个维度,并构建了 6 种泛化划分策略。

对 20 种特征-模型组合的系统评测揭示了一个统一模式:所有方法均能恢复扰动的转录组方向(Pearson delta 0.59-0.67),但没有任何方法能可靠区分不同变异的特异性响应(PDS 0.43-0.49,处于或低于随机基线 0.50)。我们将这一"方向-排序解离"现象追溯到测量分辨率地板:split-half 功效分析表明,变异特异性排序是否可达由变异-野生型效应量、变异内采样噪声与每变异有效细胞数三者的联合区制决定,当前许多设置即便在名义上较高的深度下仍落在可排序窗口附近或之下。将分析扩展至领域最常用的基因水平 Perturb-seq 图谱(Replogle 2022),我们发现在原始测序深度下、以直接 split-half 能量距离衡量,约半数(55%)的扰动不可排序,表明检测地板是分布距离度量在有限深度下的内禀属性,而非变异水平数据的特例。

AllelePerturb 为变异水平扰动预测建立了首个标准化基准,发现当前评价分辨率受限于测量而非模型,并为未来实验设计提供了每变异所需细胞数的定量处方。数据、代码和评价协议开源于 https://github.com/Boom5426/AllelePerturb。

---

## 引言

单细胞扰动预测是计算生物学最活跃的前沿之一。给定某种扰动条件,预测细胞群体的转录组响应,对药物筛选、功能基因组学和精准医学具有直接价值。随着 Perturb-seq 等高通量技术的成熟,scGen^1、CPA^2、GEARS^3、scGPT^4 等方法在基因水平扰动预测上取得了显著进展。近期,scPerturb^5 和 PerturbBench^6 等基准的建立进一步推动了方法间的标准化比较。

然而,现有方法和基准几乎全部在基因粒度上定义扰动:预测"敲低 TP53 后"或"施加药物 X 后"的细胞状态变化。这一定义隐含了一个关键简化,即同一基因的所有扰动形式产生相同的转录组后果。在基因敲除实验中这一假设大体成立,但在临床遗传学和自然变异的语境下,它遮蔽了疾病机制中最关键的一层信息。

以 TP53 为例,该基因已知有数百种错义突变,功能后果差异巨大。R175H 通过蛋白折叠丧失导致完全功能缺失,R273C 保留折叠但丧失 DNA 结合,R248Q 同时具有显性负效应和获得性功能^7。这些不同变异激活不同的转录程序,影响不同的药物敏感性和临床预后。类似地,在 KRAS 中,G12D 和 G12V 虽然都是 12 号密码子的激活突变,但对下游 MAPK 通路的激活强度存在可测量的差异^8。因此,将扰动预测的分辨率从基因水平推进到蛋白编码变异水平,是连接基因型与细胞表型的必要一步。

实现这一推进面临两个根本性问题。第一,缺乏系统性的多基因基准数据集。现有变异水平的单细胞扰动数据分散在不同实验室和技术平台上,缺少统一的处理流程和质控标准。近期的 PerturbNet^9 将变异水平扰动引入条件生成建模,但仅覆盖 GATA1 单个基因,且未提供标准化的多基因评价框架。第二,缺乏针对变异分辨率的评价协议。现有评价指标(如 MSE、Pearson 相关)衡量的是预测值与真实值之间的全局吻合度,但没有测试模型是否能正确区分同一基因内不同变异的特异性响应。这恰恰是变异水平预测的核心要求。

为解决这两个问题,本文提出 AllelePerturb,一个面向蛋白编码变异水平的单细胞扰动预测基准与评价框架。AllelePerturb-Bench 整合了 4 个基因(TP53、KRAS、GATA1、JAK1)的 472 个蛋白编码变异和 321,043 个单细胞,覆盖 Perturb-seq、碱基编辑和 scSNV-seq 三种实验技术,涵盖从 189 到 1,154 氨基酸的蛋白质长度范围。AllelePerturb-Eval 定义了包含 10 个指标的评价协议,覆盖排序准确性(PDS)、方向恢复(Pearson delta)和差异表达保真度(DE overlap)三个维度。为测试泛化边界,我们设计了 6 种划分策略:随机划分、位置外推、机制外推、跨基因转移、低深度和兼容性划分。

在 AllelePerturb 上,我们系统评测了 20 种特征-模型组合,发现了一个一致且此前未被识别的模式:所有方法均能恢复扰动的转录组方向,但不能区分不同变异的特异性响应。我们进一步证明这一"方向-排序解离"源于测量分辨率地板而非模型能力不足,并将这一发现从经验观察提升为可量化的计量学命题。这些结果重新定义了变异水平扰动预测的评价标准,并为下一代实验设计提供了每变异所需细胞数的定量参考。

---

## 结果

### AllelePerturb 建立变异水平扰动预测的基准与评价框架

现有单细胞扰动预测基准在基因粒度上定义扰动,将同一基因的所有变异归为一个条件。然而,同一基因内不同蛋白编码变异可能产生截然不同的转录组后果。以 TP53 为例,在 UMAP 空间中,基因水平的表示将所有变异折叠为单一方向(图 1a 上),而变异水平的解析揭示 R175H、R273C、R248Q 等突变沿不同轨迹偏离野生型(图 1a 下),反映了折叠丧失、DNA 结合丧失和获得性功能等不同突变机制对转录调控的差异化影响。

AllelePerturb-Bench 整合了四个具有不同生物学特征的基因(图 1b):TP53(393 aa, 98 个变异)和 KRAS(189 aa, 93 个变异)来自 A549 细胞系的 Perturb-seq 实验^10,GATA1(413 aa, 255 个变异)通过造血干祖细胞的碱基编辑实现^9,JAK1(1,154 aa, 25 个变异)来自 HT-29 细胞系的 scSNV-seq^11。四个基因共计 472 个蛋白编码变异和 321,043 个单细胞,覆盖了从 189 到 1,154 氨基酸的蛋白长度范围和从信号噪声比极低(TP53, KRAS)到较高(JAK1)的生物学信号梯度。各基因间的每变异细胞深度差异显著(图 1c):TP53 和 KRAS 的中位数分别为 929 和 658 个细胞/变异,GATA1 为 355 个,JAK1 仅 104 个,我们以 200 个细胞/变异作为参考线标注。

为编码变异的结构和生物物理学特性,我们定义了 6 维特征向量 θ(图 1d):疏水性变化、侧链体积变化、电荷变化、折叠核心位置、功能开关残基和致病性热点标记。θ 在蛋白质序列空间中编码了超越基因身份的变异性质,使不同基因的变异可在统一的特征空间中比较;在 PCA 投影中(图 1d 右),四个基因的变异呈现清晰的聚类结构,同时基因间存在重叠区域,表明 θ 捕获了跨基因的共享生物物理学信息。AllelePerturb-Eval 将评测分解为两个核心维度(图 1e):Pearson delta 衡量预测方向与真实方向的全局吻合度,PDS(perturbation discrimination score)衡量预测响应谱能否正确识别目标变异。为测试泛化边界,我们设计了随机划分、位置外推、机制外推、跨基因、低深度和兼容性 6 种划分策略(图 1e 底部)。

### 所有方法恢复扰动方向但不能区分变异身份

我们对 20 种方法(6 种模型头 × 3 种特征空间 + Gene mean 和 WT null 基线)进行了系统评测,结果揭示了一个清晰且一致的模式(图 2)。在变异排序维度上,所有 19 种非基线方法的 PDS(cosine)均落在 0.43 至 0.49 范围内,处于或低于随机基线 WT null 的 0.50(图 2a 左),最佳方法(KNN-esm+θ, 0.491)与最差方法(Gene mean, 0.437)的差距仅 0.054,PDS-L1 和 PDS-L2 给出一致结论。以 held-out 变异为重采样单位的 95% bootstrap 置信区间显示,没有任何方法的 PDS 显著超过 0.50 chance 水平:19 种方法中 16 种的区间完全落在 chance 之下,其余 3 种与 chance 重叠(图 2a 左须线)。与之形成鲜明对比的是,在方向恢复维度上,所有方法的 Pearson delta 均显著高于零基线,范围为 0.59 至 0.67(图 2a 右),最佳方法 Lasso-esm 达到 0.669,且全部 19 种方法的 Pearson delta 置信区间均完全高于零。将两个维度置于同一散点图中(图 2b),全部方法均落在左上象限(方向恢复但排序失败),理想区域空无一人,且这一模式在所有模型类型和特征空间中完全一致,表明它是当前评价分辨率下的系统性现象。

差异表达保真度提供了对这一解离的生物学注释(图 2c):direction agreement 78%,DE-LFC Spearman 0.46,DE overlap 仅 28%。这一梯度(方向 78% > LFC 排序 46% > 基因集 28%)表明模型学到了扰动的"程序级"方向,但没有捕获使不同变异可区分的精细信号。

按基因分解(图 2d),解离幅度呈现有规律的梯度:TP53 最强(Δ = 0.31),其次是 KRAS(0.18)和 GATA1(0.15),JAK1 最弱(0.02)。这一梯度与信噪比一致:JAK1 逐变异信号最强(split-half D_self/D_null = 0.14),解离最弱;TP53/KRAS 信号接近纯噪声(D_self/D_null 分别为 0.96 和 1.01),解离最强。在单个变异尺度上(图 2e),以 TP53 Ridge-esm 为例,每个 held-out 变异的 Pearson delta 稳定在 0.74 至 0.82,但 PDS 从 0.00 到 0.88 剧烈波动,8 个代表性变异中有 4 个低于 chance。这一规律一致地提示:排序失败并非源于模型能力不足,而可能与测量分辨率有关。

### 方向-排序解离源于测量分辨率地板

上述解离的一个自然假说是:变异间的转录组差异确实存在,但在当前每变异细胞数下无法被分布距离度量稳定检出。为验证这一假说,我们设计了 split-half 可复现性检验。对每个变异,将其细胞随机等分为两半,计算两个半集之间的能量距离(D_self)以及野生型与变异半集之间的能量距离(D_null)。若 D_self/D_null 接近 1,则意味着变异内部的采样波动与变异-野生型之间的真实信号相当,排序在统计上不可达(图 3a)。

四个基因的 D_self/D_null 比值揭示了一个清晰的梯度(图 3b):TP53(0.96,95% CI 0.95–0.98)、KRAS(1.00,0.99–1.02)和 GATA1(0.88,0.86–0.90)均接近 1,表明在当前细胞数下,这三个基因的变异特异性信号在很大程度上被采样噪声淹没;JAK1(0.21,0.12–0.33)显著低于 1 且置信区间排除了噪声地板,对应其更强的逐变异效应和更深的有效检测余量。GATA1 略低于噪声线,并带有一条个体可排序变异的长尾,与其中等解离幅度一致。这一梯度与图 2d 的解离幅度吻合:信号被噪声淹没的基因,解离最强。

因此,一个变异是否可排序并非仅由细胞数决定,而是由效应量、变异内采样噪声与每变异有效细胞数三者的联合区制决定。将每个基因按其变异-野生型中位效应量对 split-half 比值作图使这一点显式化(图 3c):JAK1 的中位效应量比其余基因高一个数量级,因而仅需每半 50 个细胞即可完全可检测;而 TP53 和 KRAS 即便在名义深度高得多(中位 929 和 1,000 原始细胞/变异)时仍停留在噪声地板。在 n = 50 至 300 细胞/半的范围内对 split-half 检测率作子采样验证了这一排序(图 3e):JAK1 的检测率在整个范围内饱和于 100%,而 TP53、KRAS、GATA1 仅随深度小幅上升并在 75% 以下趋于平台。只有当效应窗口足够宽、增加的深度能够分辨它时,增加细胞才有帮助;对低效应基因而言,是窗口本身而非细胞数构成了约束。这解释了为什么 JAK1 的 PDS 接近理想区而其余三个基因均低于 chance。

这一检测地板是否仅为变异水平扰动数据的特例?为回答这一问题,我们将 split-half 可排序性分析扩展至领域最常用的三个基因水平 Perturb-seq 图谱(图 3d),全程使用一个统一定义:原始最大深度、全部扰动、直接 split-half 能量距离。四个变异水平基因按此定义有 98% 至 100% 不可排序,JAK1 作为唯一的阳性对照为 10%。在基因水平图谱中,Replogle 2022(K562, 1,832 个扰动,中位 132 细胞/扰动)有 55.3%(95% CI [52.8%, 57.5%])的扰动不可排序,Norman 2019(236 个扰动)为 3.4%,Adamson 2016(96 个扰动)为 14.6%,该梯度与各图谱的逐扰动深度一致。结果在显著性过滤和匹配深度的分析变体下定性稳健,尽管绝对比例随聚合选择而变(Extended Data)。这一发现为近期文献中报道的基因水平扰动预测方法间差距极小的现象提供了计量学解释:当大比例扰动落在评价指标的分辨率之下时,方法间的排名本身就不稳定。


### 方向-排序解离在不同划分策略、距离函数和特征空间中一致成立

上述发现可能受特定划分策略、距离度量或特征选择的影响。为排除这些可能,我们在 5 种泛化划分、3 种 PDS 距离函数和 3 种特征空间上进行了全网格评测(图 4)。

PDS 和 Pearson delta 在不同划分下表现出截然不同的稳定性(图 4a, b)。所有方法的 Pearson delta 在随机、位置外推和兼容性划分下保持 0.55 至 0.75 的正值,仅在机制外推和低深度划分下有所下降。PDS 则在所有划分下均处于 0.29 至 0.54 的 chance 带内,低深度划分尤为严重(最佳方法 PDS = 0.36)。排序失败不受距离函数选择的影响(图 4c):PDS-cosine(中位 0.464)、PDS-L1(0.476)和 PDS-L2(0.474)在 18 种方法上的分布完全一致,全部低于 chance。

将分析从方法维度转向基因-划分组合(图 4d),17 个基因 × 划分组合中有 15 个呈现解离模式:除两个机制外推组合(JAK1 与 KRAS split 3,其 held-out 集较小、两项指标接近)外,Pearson delta 均高于对应的 PDS。GATA1 Low-N 划分的 PDS 降至 0.19,验证了图 3 的检测极限预言。三种特征空间的比较进一步排除了表征能力的解释(图 4e):ESM-1v(1,280 维)的 Pearson delta(0.64)与 6 维 θ(0.62)相当,两者拼接后(0.63)无增益;三者的 PDS 均为 0.45。从 6 维手工特征到 1,280 维蛋白语言模型嵌入,方向恢复已饱和,而排序瓶颈不在特征维度。

现有单细胞扰动预测工具同样无法适配变异水平任务(图 4f)。CPA 的 dosage 接口仅接受标量输入,无法传入多维 θ 向量;CellFlow 和 scDFM 缺乏对训练集外变异的预测接口;GEARS 以基因为图节点,在架构层面对变异不可见;STATE 需要对每个基因独立训练且在 GATA1 规模上内存溢出;唯一成功运行的 Biolord 将 θ 作为有序属性输入,但 PDS 仍落在 chance 附近(0.47 至 0.51)。这些失败模式指向一个共同原因:现有方法假设扰动身份在基因或药物水平定义,缺乏显式的变异条件化生成接口。

### 可排序性可先验预测,但绝对校准依赖数据集

检测极限分析不仅解释了当前评价的系统性失败,还提供了一个实用工具:在建模之前预判哪些变异-数据集组合处于可评价窗口内。我们构建了一个 leave-one-dataset-out(LODO)可排序性预测器,以效应量、每变异细胞数和基因空间维度等先验可获取的特征预测每个扰动是否可排序(图 5a, b)。

该预测器在 7 个数据集的 LODO 交叉验证中展现出优异的判别能力(图 5a):核心 4 个数据集的平均 AUROC 为 0.974,全部 7 个数据集为 0.965。JAK1(AUROC = 0.999)和 Replogle(0.981)的判别几乎完美,即使是最难的 TP53(0.943)也远超随机。这意味着在开展变异水平扰动预测之前,仅凭 split-half 效应量估计即可可靠地筛选出哪些变异处于可评价窗口内。然而,该预测器的绝对校准能力较弱(图 5b):预测的不可排序比例与观测值之间的 R² 仅为 0.11,平均绝对误差 26 个百分点。GATA1 的偏差尤为严重(预测 23%, 观测 96%)。这一差距源于数据集间效应量分布的系统性差异,提示该预测器适用于分诊(triage)而非绝对的跨数据集校准。

基于检测极限分析,所需深度最好理解为在特定效应量与噪声区制下的条件量,而非单一的通用每变异细胞数阈值(图 5c)。JAK1 等强效应基因(变异-野生型中位效应量约 13)仅需每半 50 个细胞即可突破检测地板;而 TP53(效应量约 1.4)和 KRAS(约 1.1)即便在中位原始深度 929 和 1,000 细胞/变异下仍停留在噪声地板(D_self/D_null 0.96 和 1.01)。因此,每变异细胞数目标应读作特定效应量、变异间差异与采样噪声结构下的一个功效估计,而非固定标准:TP53 和 KRAS 名义深度高,但等位特异性效应窗口窄,故仍接近排序地板。这将瓶颈从"细胞数不足"重新框定为"效应量-噪声窗口不足"。

将同一套统一定义应用于领域最常用的基因水平 Perturb-seq 图谱(图 5d),在原始深度下 Replogle 有 55% 的扰动不可排序,Adamson 为 15%,Norman 为 3%,该排序与各图谱的逐扰动深度一致。当归一化到相同的浅深度(每半 n = 50)后,三者的比例相互靠拢(Replogle 55%、Adamson 39%、Norman 11%):深采样图谱基本不变,而较浅的两者恶化,证实检测地板是分布距离度量在有限深度下的内禀属性,而非任一特定数据集的特性。这一结果为近期文献中报道的基因水平扰动预测方法间差距极小的现象提供了计量学解释。

综合以上发现,我们提出了一个功效感知的变异水平扰动评价工作流(图 5e):先通过小规模先导实验(10 至 50 个变异,每变异约 100 个细胞)估计效应量和 split-half 地板;再利用可排序性预测器筛选可评价的变异;据此设计实验,使效应量-噪声窗口(而不仅是原始细胞数)能够突破检测地板;最后仅在信号超过该地板的变异上评价预测模型。这一工作流将"所有变异一视同仁地评价"转变为"功效感知的分层评价",避免在测量分辨率不足的条件下产生误导性的方法排名。


## 讨论

本文通过 AllelePerturb 基准的系统评测,揭示了单细胞扰动预测领域一个此前未被识别的评价分辨率瓶颈:所有测试方法均能恢复扰动的转录组方向,但不能区分同一基因内不同变异的特异性响应。这一"方向-排序解离"并非源于模型能力不足,而是分布距离度量在当前数据深度下的内禀测量极限。

这一发现为近期文献中的多个观察提供了统一解释。Ahlmann-Eltze 等在 Nature Methods 2025 中报告,基因水平扰动预测的深度学习方法尚未系统性超越简单线性基线,顶方法间的性能差距极小。我们的检测极限分析表明,当评价指标本身的分辨率不足以区分扰动时,方法间的排名必然不稳定,与方法复杂度无关。类似地,ADAPERT 报告的 scVI PDS 0.502(标准差 0.290)和 Chreode 排行榜上方法间仅 0.014 的 Pearson 相关差距,均与"评价处于检测地板附近"的解释一致。我们的贡献不是重复这些观察,而是将其追溯到一个可量化的计量学根源,并提供了诊断工具(split-half 可排序性检验)和实验设计处方(每变异所需细胞数)。

方向-排序解离的生物学含义值得强调。方向恢复成功意味着不同变异确实共享一个"程序级"的转录组响应方向:TP53 错义突变整体偏向功能丧失程序,KRAS 激活突变偏向 MAPK 程序。这一共享方向可被简单的线性模型从 6 维生物物理特征中捕获。但排序失败意味着使不同变异(如 R175H vs R273C)可区分的那部分精细信号,在当前每变异数百个细胞的规模下不可稳定测量。这不是说变异间没有转录组差异(文献中有充分的生化证据),而是说当前单细胞扰动实验的统计功效不足以在分布距离框架下检出这些差异。

将分析扩展至基因水平 Perturb-seq 图谱(Replogle、Norman、Adamson)后,我们发现检测地板是分布距离度量的内禀属性,而非变异水平数据的特例。在 Replogle 图谱的原始深度下,以直接 split-half 能量距离衡量约半数(55%)的扰动不可排序,且三个图谱的不可排序比例与其逐扰动深度一致,在匹配的浅采样下相互靠拢。这进一步确认了地板的计量学本质。这意味着领域中使用 E-distance 或类似分布距离作为评价指标的所有基准,都可能在低深度区间产生不稳定的方法排名。

本研究存在若干局限。第一,AllelePerturb 目前仅覆盖 4 个基因,且变异类型以错义突变为主,对移码、无义和剪接变异的覆盖有限。第二,6 维 θ 特征是手工设计的,虽然在方向恢复上与 1,280 维 ESM 嵌入相当,但可能遗漏了某些变异属性。第三,split-half 检验依赖于细胞数的对半分割,在极低深度(少于 50 个细胞)时自身的统计功效受限。第四,LODO 可排序性预测器的判别能力优异(AUROC 0.97)但绝对校准较差(R² = 0.11),限制了其作为跨数据集定量工具的应用。第五,我们未能成功运行 CPA、CellFlow、scDFM 等外部方法的变异水平评测,这些方法的失败源于接口不兼容而非方法本身的缺陷。

展望未来,突破检测极限需要在两个方向上推进。在实验设计端,下一代变异水平扰动实验(如 Tahoe-100M 规模的碱基编辑图谱)有望提供远超当前深度的每变异细胞数,使变异特异性排序进入可评价窗口。在方法设计端,变异水平扰动预测需要显式的变异条件化生成接口,而非将变异身份映射到现有基因/药物水平的扰动框架上。AllelePerturb-Eval 的功效感知评价工作流为这两个方向提供了标准化的评价基础设施。

---

## Methods

### Benchmark construction

**Data sources.** AllelePerturb-Bench integrates single-cell perturbation data from four genes across three experimental technologies. TP53 (98 variants, 83,422 cells) and KRAS (93 variants, 83,556 cells) were obtained from Perturb-seq experiments in A549 cells (GSE161824, Ursu et al. 2022). GATA1 (255 variants, 149,161 cells) was obtained from base-editing screens in human hematopoietic stem and progenitor cells (HSPCs), processed from PerturbNet (Yu and Welch 2025, HuggingFace: cyclopeta/PerturbNet_reproduce). JAK1 (25 variants, 4,904 cells) was obtained from scSNV-seq in HT-29 cells treated with IFN-gamma (PRJEB48915, Cooper et al. 2024).

**Preprocessing.** For TP53 and KRAS, we loaded the processed count matrices, barcodes, and variant-to-cell mappings from GEO. Cells were assigned to variants using the published variant-to-cell CSV files. For GATA1, we used the PerturbNet-preprocessed AnnData object containing log-normalized expression over 2,477 highly variable genes (HVGs) with predefined holdout splits. For JAK1, we loaded the SingleCellExperiment object in R, filtered to homozygous genotyped cells with clean consequence annotations (missense, splice-site, stop-gained, synonymous), selected the top 2,000 HVGs using scran::modelGeneVar(), and exported log-normalized counts. All datasets were z-scored per gene across cells within each dataset. Wild-type (WT) cells served as the shared reference within each dataset.

**Variant feature vector (theta).** Each protein-coding variant was encoded as a 6-dimensional biophysical feature vector theta = (delta_hydrophobicity, delta_sidechain_volume, delta_charge, fold_core_location, functional_switch_residue, hotspot_pathogenic_marker). The first three features were computed as the difference between mutant and wild-type amino acid properties using the Kyte-Doolittle hydrophobicity scale, amino acid molecular weights, and formal charge at pH 7.0, respectively, then z-scored across all 472 variants. fold_core_location was set to 1.0 for residues within known structural core domains (TP53 DNA-binding domain 94-293, KRAS G-domain 1-169, GATA1 zinc fingers ZF1 204-228 and ZF2 258-282, JAK1 kinase domain 866-1154), 0.5 for interaction regions, and 0.0 elsewhere. functional_switch_residue was set to 1.0 for known catalytic or coordination residues (e.g., KRAS G12/G13/Q61, GATA1 zinc-finger cysteines). hotspot_pathogenic_marker was defined using external prior annotations only, with no reference to the single-cell outcome: canonical COSMIC/IARC missense hotspot codons for TP53 (175, 176, 179, 220, 238, 245, 248, 249, 273, 282, 285) and KRAS (12, 13, 59, 61, 117, 146), zinc-finger zinc-coordinating cysteines for GATA1 (204, 207, 225, 228, 258, 261, 279, 282) together with ClinVar pathogenic/likely-pathogenic annotations, and functional-module membership (JH2 pseudokinase 583-855, JH1 kinase 875-1153) for JAK1. An earlier version of this feature used the top 20% of measured effect size, which introduces outcome leakage; we removed that term and re-ran the full evaluation grid. Ranking performance was unchanged (theta-only PDS shifted by at most 0.02, and all methods remained at or below the chance baseline), confirming that the direction-ranking dissociation does not depend on any outcome-derived feature.

**ESM embeddings.** For ESM-based feature representations, we extracted per-variant embeddings from ESM-1v (esm1v_t33_650M_UR90S_1) by encoding the full mutant protein sequence and taking the mean of the last hidden layer over all residue positions, yielding a 1,280-dimensional vector, which was used directly as the ESM feature without dimensionality reduction.

### Evaluation protocol (AllelePerturb-Eval)

**Pseudobulk perturbation profiles.** For each variant v, we computed the pseudobulk perturbation effect as delta_v = mean(X_v) - mean(X_WT), where X_v and X_WT are the z-scored expression matrices of variant and wild-type cells, respectively.

**Perturbation discrimination score (PDS).** PDS measures whether a method's predicted perturbation profile for variant v is closer to the true profile of v than to those of other variants. For each held-out variant v, we computed the distance d(pred_v, true_v) and ranked it among {d(pred_v, true_w) : w in all test variants}. PDS = 1 - (rank - 1) / (n_test - 1), where rank 1 (correct match is nearest) gives PDS = 1 and random matching gives E[PDS] = 0.5. We computed PDS using three distance functions: cosine distance (PDS_cos), L1 distance (PDS_L1), and L2 distance (PDS_L2).

**Pearson delta.** The Pearson correlation between the predicted pseudobulk delta and the true pseudobulk delta, averaged across held-out variants.

**Delta cosine.** The cosine similarity between predicted and true pseudobulk deltas.

**Pearson delta top-20.** Pearson delta restricted to the 20 genes with highest variance across all variants, testing whether signal concentrates in a small gene subset.

**DE overlap.** For each variant, we identified the top 50 differentially expressed genes (by absolute delta) in both predicted and true profiles. DE overlap = |DE_pred intersect DE_true| / 50.

**DE-LFC Spearman.** Spearman rank correlation of log-fold-change values restricted to the union of predicted and true top-50 DE gene sets.

**Direction agreement.** Fraction of the top-50 DE genes (by true absolute delta) for which the predicted and true deltas have the same sign.

**MAE.** Mean absolute error between predicted and true pseudobulk deltas, averaged across genes and variants.

### Generalization splits

We defined six split strategies to test different aspects of model generalization:

1. **Random (Split 1):** 35% of variants held out uniformly at random per gene, using predefined holdout assignments where available (GATA1 holdout1).
2. **OOD-Position (Split 2):** Variants in the C-terminal half of each protein held out for testing.
3. **OOD-Mechanism (Split 3):** Variants annotated as hotspot or functional-switch residues held out for testing.
4. **Cross-Gene (Split 4):** All GATA1 variants held out; train on TP53+KRAS+JAK1 only. (Not included in the main grid due to gene-space incompatibility.)
5. **Low-N (Split 5):** Variants with fewer than 200 cells held out for testing.
6. **PerturbNet-Compatible (Split 6):** Held-out sets matching the PerturbNet evaluation protocol for GATA1.

### Prediction methods

We evaluated 20 feature-model combinations constructed from 6 model heads and 3 feature spaces, plus 2 baselines:

**Model heads.** (1) Ridge regression (alpha=1.0), (2) Lasso regression (alpha=0.01), (3) Random Forest (100 trees), (4) Gradient Boosted Trees (GBoost, 100 estimators), (5) K-nearest neighbors (KNN, k=5), (6) Multi-layer perceptron (MLP, 2 hidden layers of 64 units, ReLU, 500 epochs). All models were trained to map variant features to pseudobulk delta profiles.

**Feature spaces.** (1) theta: 6-dimensional biophysical features, (2) ESM: 1,280-dimensional ESM-1v embeddings, (3) ESM+theta: concatenation of theta and ESM features (1,286 dimensions).

**Baselines.** (1) WT null: predicted delta = 0 for all variants (no perturbation effect). (2) Gene mean: predicted delta = mean delta across all training variants of the same gene.

**External methods.** We attempted to evaluate six external perturbation prediction methods: CPA, CellFlow, STATE, GEARS, scDFM, and Biolord. CPA's dosage interface accepts only scalar inputs and produced NaN losses when provided with multi-dimensional theta vectors. CellFlow and scDFM lack prediction interfaces for variants not seen during training. GEARS uses a gene-keyed knowledge graph that is allele-blind by design. STATE required per-dataset retraining and exceeded GPU memory on GATA1 (149K cells). Biolord was the only external method that ran successfully, using theta dimensions as ordered attributes; its PDS ranged from 0.47 to 0.51 (chance level).

### Split-half rankability analysis

To determine whether allele-specific ranking is measurement-limited, we performed split-half reproducibility analysis. For each variant with at least 50 cells, we randomly split its cells into two equal halves (A and B) and computed:

- D_self = E-distance(A, B): the within-variant replicate distance
- D_null = E-distance(A, WT_sample): the variant-to-wild-type distance

where E-distance(X, Y) = 2 * mean(d_XY) - mean(d_XX) - mean(d_YY) with Euclidean distances in PCA-50 space, subsampling n=min(n_half, 150) cells per group. A variant was classified as "rankable" if its signal S = D_null - D_self exceeded the 95% bootstrap confidence interval width W of D_self (i.e., S > W), and "un-rankable" otherwise. This procedure was repeated with 50 bootstrap resamples to estimate confidence intervals.

The analysis was applied to all four AllelePerturb genes and extended to three gene-level Perturb-seq benchmarks: Replogle 2022 (K562, 2,003 perturbations), Norman 2019 (K562, 236 perturbations), and Adamson 2016 (K562, 103 perturbations), downloaded from the scPerturb data portal.

### LODO rankability predictor

We trained a logistic regression classifier to predict per-perturbation rankability using leave-one-dataset-out (LODO) cross-validation across all 7 datasets. Features included: log(cells per perturbation), pseudobulk effect size (L2 norm of delta), PCA explained variance ratio, and gene-space dimensionality. The classifier was evaluated by AUROC (discrimination) and by comparing predicted versus observed un-rankable fractions per held-out dataset (calibration).

### Power curve analysis

To estimate the number of cells per variant required for reliable detection, we subsampled cells at n = {50, 100, 150, 300} per variant and recomputed the split-half rankability test at each depth. The fraction of variants classified as "detectable" (rankable) was plotted as a function of subsample size for each gene. Detection thresholds (cells at 50% detection rate) were estimated by sigmoid curve fitting where sufficient data points were available.

### Statistics

Confidence intervals for PDS and Pearson delta were estimated by nonparametric bootstrap resampling with the held-out variant as the resampling unit (2,000 resamples, percentile method, 95% level); a method was called significantly below chance when the upper bound of its PDS interval fell below 0.50, and significantly positive on direction when the lower bound of its Pearson-delta interval exceeded 0. The same bootstrap was applied to the per-gene D_self/D_null ratio (resampling over variants) to obtain the intervals reported in Fig. 3b. Rankability confidence intervals for the external Perturb-seq atlases (Fig. 3d) were computed by bootstrapping over perturbations. All bootstraps used a fixed random seed for reproducibility. No adjustment for multiple comparisons was applied, as the conclusions rest on the direction and consistency of effects across all methods rather than on any single significant test.

### Software and reproducibility

The evaluation grid was run in Python 3.11 using NumPy 2.4.6, pandas 2.3.3, scikit-learn 1.7.2, SciPy 1.14.1, and matplotlib 3.10.8. Single-cell data loading and preprocessing used scanpy 1.10 (Python) and scran with SingleCellExperiment (R, Bioconductor). ESM embeddings were extracted with fair-esm 2.0.0 and PyTorch 2.10. Figures were rendered with matplotlib 3.10. All code and evaluation protocols are available at https://github.com/Boom5426/AllelePerturb.

