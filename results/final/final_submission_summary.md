# Final Submission Summary

最终投稿目录：

`D:\AAUDE\paper\paper4\revision_v2\fpy\submit_version`

## 有效文件

1. `manuscript_260309.docx`
   - 标题：`Asymmetric impacts and divergent futures of drought and waterlogging on global crop yields`
   - 结构：227 paragraphs，189 non-empty paragraphs，0 tables，1 section。
   - 核心主题：通过 CWatM 内嵌反事实框架，归因 drought 与 waterlogging 对全球主要作物产量损失的非对称影响。

2. `SI_20260228.docx`
   - 结构：150 paragraphs，121 non-empty paragraphs，2 tables，1 section。
   - 核心内容：Supplementary Methods，详细说明 root-zone soil moisture stress factor、aeration stress、水分-产量响应、counterfactual framework、验证和参数化。

3. `comments_response.docx`
   - 结构：179 paragraphs，175 non-empty paragraphs。
   - 核心内容：回应三位审稿人，重点解释从旧版 “coupled water and heat stress / extreme climate” 收缩为 “drought and waterlogging attribution” 的原因。

4. `cover letter.docx`
   - 结构：11 paragraphs，10 non-empty paragraphs。
   - 核心内容：投稿 Nature Communications，突出 process-based attribution、drought-waterlogging asymmetry、国家粮食系统风险。

## 最终稿与旧稿的关键区别

最终稿不再以 `coupled extreme water and heat stress` 为主线，而是明确聚焦：

- drought vs waterlogging
- abnormal water stress
- embedded counterfactual yield experiment
- process-consistent attribution
- future divergence under SSP1-2.6 and SSP5-8.5
- national food-system risk

这意味着后续修改时，不能再沿用旧版中关于 heat stress、extreme precipitation 或 Dr index 的强表述，除非用户明确要求恢复并补充证据。

## 当前 Methods 核心定义

- Drought：SPEI < -1 且 root-zone soil moisture < crop-specific critical threshold。
- Waterlogging：SPEI > 1 且 root-zone soil moisture > field capacity。
- SPEI 只用于识别异常气候状态，不直接进入 crop water-stress 或 yield calculation。
- Yield loss 通过 actual/stressed simulation 与 counterfactual normal water-stress simulation 的相对差异定义。
- Irrigated/rainfed 需要区分：模拟系统与空间面积权重，不应被写成灌溉实践因果效应，除非代码确实支持该解释。

