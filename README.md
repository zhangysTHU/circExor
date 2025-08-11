# circExor

#### 介绍
目前已有大量分子生物学的证据指出，RNA本身的序列信息便涵盖了大量定位信息，多种RNA的核苷酸序列信息可以作为顺式调节原件，调控RNA的胞外分选或决定RNA的细胞内亚定位。这为使用机器学习的手段通过序列预测非编码RNA定位提供了可行性基础。目前已有的研究诸如使用深度学习模型预测和解释各类sRNA向胞外的分选，使用梯度提升树解析long RNA在细胞中的亚定位等。环状RNA（circRNA）作为一种重要的非编码RNA分子，已被发现参与多种生物学过程，包括基因表达调控和细胞内信号传递 。近年来，研究发现circRNA不仅存在于细胞内，还可以通过外泌体（EV）释放到细胞外，成为细胞间信息传递的媒介。我们注意到，目前的纯计算定位研究主要集中于long RNA、microRNA领域，学习目标主要为解析非编码RNA的亚细胞定位，对RNA的细胞外定位/分选涉猎较少，特别是circRNA，目前尚没有针对其序列信息的机器学习研究。我们从RNAlocate数据库中筛选和导出了大量带有cell-EV定位标签的circRNA数据，并根据circRbase等数据库构建了完整的circRNA序列-定位数据集，通过机器学习的手段在自建数据集上训练了名为circExor的胞内-EV定位预测工具，首次为预测circRNA胞外定位预测提供了低成本纯计算工具和首创数据集，同时进一步验证了核酸序列作为RNA外泌的顺式调节原件。

#### 研究方法
1.**circRNA location dataset construction**：
    RNALocate 是一个可通过网络访问的数据库（http://www.rna-society.org/rnalocate/），提供了多个通过实验证据得到或通过已发表的机器学习模型得到的RNA定位条目 。我们从数据库中导出仅基于实验证据的，circRNA类别的所有RNA条目，包括来源于人类和小鼠的134184条circRNA数据。RNAlocate数据库提供了每条circRNA的RNA_Symbol，通过已发表的实验结果推测出的RNA定位位置（Subcellular_Localization），和数据库作者提供的定位可靠性分数（RNALocate_Score）。根据RNA_Symbol，我们从circbase 和exoRbase 数据库的fasta文件导出相关circRNA的cDNA序列，包含人类来源circRNA，小鼠来源circRNA，外泌体circRNA等。无法从以上数据库检索到序列的circRNA定位条目舍去。合并上述文件后我们得到包括RNA_Symbol，Subcellular_Localization，RNALocate_Score，cDNA sequence的数据框。
在随后的数据清理中，为保证定位信息的可靠性，我们清除了所有RNALocate_Score≤0.5的条目，并删除了所有序列中含有N的条目。出于机器学习工程上的考量，我们限制circRNA序列长度在0 nt到11238 nt之间，并删除了若干RNA数量极为稀少的定位位置或特殊定位的circRNA条目。至此共有737 个细胞核（Nucleus）定位条目，168个细胞质（Cytosol）定位条目，和385 个EV（Extracellular vesicle）定位 cricRNA，共1290个circRNA定位条目，其中细胞核和细胞质定位合并为胞内定位905条。
2.**Calculation of circRNA k-mer frequency calculation and machine learning**： 
    在进行k-mer编码前对两类样本进行随机下采样，并以胞内定位circRNA作为负类，EV定位circRNA作为正类。参考目前已有研究中对RNA序列的编码方式，我们选择k-mer frequency量化每条circ RNA的序列特征。我们从k = 3,4,5出发计算，将每条cDNA序列转化为1344种短ATCG核苷酸组合，随后分别计算所有样本中每个k-mer的出现频率，将每个cDNA序列计算成长为1344的一维向量，向量中每个数值为一个k-mer的出现频率，为0到1之间的浮点数。1290个circRNA定位条目转化为1290*1344的矩阵作为机器学习输入。
    我们选择了逻辑回归（LR），支持向量机（SVM），随机森林（RF），lightGBM作为候选的机器学习模型，每个模型使用上文得到的相同输入矩阵。四个模型的超参列表可参考**仓库文件架构**项目中的相关代码。
    在每个模型训练过程中，我们将所有样本随机划分20%作为测试集，其余作为训练集，在通过五折交叉验证和随即参数搜索的方式于训练集上调参，使用训练集上的最优参数在测试集上验证学习效果。
3.**Identification of EV and Intracellular resident related k-mers**：
    树 SHAP值用于测量 k-mer 对 RNA 亚细胞定位的贡献。在SHAP方法中，根据最佳信用分配的博弈论 Shapley 值计算的 SHAP 值被分配给树模型中特定 RNA 的每个k-mer 。对于每个与EV和胞内相关的 k-mer，根据在不同SHAP值与k-mer frequency之间的相关性，可以判断每个k-mer对其所在circRNA定位的影响，如SHAP值与k-mer frequency高正相关的k-mer表示其倾向于提供EV定位分类，反之表示其倾向于提供胞内驻留定位分类；而平均绝对 SHAP 值则用于提供对任意一种k-mer的 RNA 定位重要性。


#### 仓库文件架构
本仓库提供了circExoer项目的全部训练评估代码，和部分训练资料。
1.  `/SHAP/`:存放对于机器学习结果进行shaply计算和绘图的脚本和结果
2.  `/_plot/`:存放用于md文档中的图片
3.  `/models/`:构建机器学习模型和训练的代码，以及保存的机器学习结果和机器学习输出
4.  `/resources/`:自locate数据库中初筛的circRNA定位信息
5.  `/sample_preprocessings/`:对/resources/中资源和人类基因组fasta文件进行预处理的脚本


#### 未来开发代办
- 美化sample_preprocessings中涉及的数据集分析统计图
- 添加/motif分析部分，使用MEME suite中的TomTom对
- 在evaluation中添加circlocnet作为新的基线对比模个候选kmer进行motif富集
- 英文版README


## License

This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.

This repository contains modified code from the [RNAlight](https://github.com/YangLab/RNAlight.git) project, which is licensed under the MIT License.

