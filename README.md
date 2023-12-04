# Link Prediction

Study bias, topological imbalance, and biased predictive models are cyclicly connected, creating a vicious cycle that hinders our ability to make accurate predictions about less studied genes, diseases, and drugs.

It has recently been shown how study bias in protein binding experiments creates topological imbalanced networks (Lucchetta et al. 2023). Contemporary research by (Bonner et al. 2022) showed that knowledge graph embedding (KGE) link prediction (LP) models, when trained on topological imbalance networks become heavily biased towards recommending nodes with high degrees (Bonner et al. 2022). We show these overly prioritized nodes are those that have been extensively studied; when these predictions are used to generate hypotheses and direct experimental studies, creating more study bias. The cycle goes on and on, creating a system of preferential attachment, where well-studied (and subsequently connected) genes get more connections while less-studied nodes receive few to no more connections.

<p align="center"><img src="https://github.com/MSBradshaw/LinkPrediction/blob/main/genome_informatics_poster.png?raw=true" width="90%"/></p>

