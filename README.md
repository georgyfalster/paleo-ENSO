This script attempts to replicate Figure 3 of the following paper:

Freund, M. B., D. C.Verdon-Kidd, K. J.Allen, and J. R.Brown. 2025. “El Niño Southern Oscillation Reconstructions During the Last Millennium.”
*Wiley Interdisciplinary Reviews: Climate Change* **16(6)**, e70036. https://doi.org/10.1002/wcc.70036.

**The purpose**: that figure is similar to an analysis I had done some years previously, whilst working on <a href="[https://www.instagram.com/p/DPNuEG2EgOi/](https://www.nature.com/articles/s41586-023-06447-0)" target="_blank">this paper</a>. But the results were quite different. I was curious as to why. 

I contacted the authors about this (initially in February 2026), and we have had a brief discussion. The authors performed the hierarchical agglomerative clustering analysis that underpins Fig. 3 in matlab (scripts not publically available), but were unable to remember precisely their data processing steps. Hierarchical agglomerative clustering is a fairly straightforward process, such that if the analysis approach that underpins Section 5.1 of Freund et al. (2025) is robust, then replicating the analysis in R (as I have done here) should give a similar result. 

In their review, Freund et al. state "*The Falster et al. (2023) Pacific Walker circulation reconstruction stands out as highly distinct from all ENSO reconstructions in the dendrogram. It forms its own branch at a much higher linkage distance, indicating low covariance similarity. This distinct clustering suggests that the Walker circulation represents a related but different aspect of tropical climate variability compared to ENSO reconstructions.*"

In the first instance, I wondered if this in fact this feature of Freund et al.'s dendrogram existed simply because the Falster et al. reconstruction targets the Walker Circulation rather than the El Niño Southern Oscillation. The indices for the two are anti-correlated by convention. If Freund et al. neglected to account for this, it could have led to an (erroneous) apparent disparity between the PWC reconstruction and the various reconstructions with an ENSO-like target. 

This does indeed seem to be the case: when all reconstructions are aligned to correlate positively, as *should* be the case for this analysis, no single reconstruction stands out as highly distinct from the others. Rather, the reconstructions simply group with others from the same study, and to a lesser degree near others using similar proxy types. 

If including all reconstructions in their 'native' direction (i.e., PWC and SOI reconstructions anticorrelated with ENSO reconstructions), the results are as expected: the PWC and SOI reconstructons cluster together. Note that this does not therefore imply that the SOI and PWC represent a different aspect of climate variability, but simply that the indices are inversely correlated by convention. Despite using identical data, an identical time period, and the same method (hierarchical agglomerative clustering), I was unable to exactly reproduce Freund et al., Figure 3. 

Finally, some details on my scripts & the figures in this repo: using `TP_recons_cluster.R` you will be able to reproduce the figure `tropicalPacific_reconstructions_1727–1900hclust_dendrogram_excludingLiu2024.png` as this uses only publicly available data. However, Freund et al. also included a reconstruction from Liu et al. (2024) which is only available on request to the authors. The results from that analysis, which are directly comparable with Freund et al.'s Fig. 3, are in `tropicalPacific_reconstructions_1727–1900hclust_dendrogram_includingLiu2024.png`. The two are extremely similar.

The version that doesn't account for the PWC/ENSO index conventions, and is most similar to Freund et al.'s Fig. 3 is `tropicalPacific_reconstructions_1727–1900hclust_dendrogram_includingLiu2024_reconsNotAligned.png`. 

In the script, I have distilled each reconstruction's metadata based on each paper's Methods section. Metadata include:
* Name of the lead author
* Publication year
* Proxy type/s (e.g., tree ring width, tree δ<sup>18</sup>O, coral δ<sup>18</sup>O etc etc)
* Reconstruction target (e.g., Niño 3.4, Nino 4, the ΔSLP index for the strength of the PWC etc etc)
* The geographical distribution of the proxy network
* Whether or not the reconstruction was nested, i.e., with sample depth that changes through time
* The method used to reconstruct the target index from the proxy data
* The interval over which the target index was averaged (e.g., a calendar year, tropical year, particular season)
* Broadly, the variables reflected by the proxies (e.g., moisture, temperature)
* Any preprocessing steps applied to the proxy data (e.g., low-pass filtering, detrending)

You can choose which of these you show on your own dendrograms. 
