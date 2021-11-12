# Model complexity affects species distribution projections under climate change

- PAPER DE REFERENCIA FOR SMDS, VERY GOOD.

## General notes

- I have NOT used this paper to support our consideration of different modelling choices, 
	- because they do not just change SDM algorithms and climate change scenarios, also vary the number of variables, the parameter complexity in each type of model, etc... We did well, but it is still not at the level of this paper. 
	- in addition, i do not fully understand their approach to estimate the importance of each factor of choice (climate scenario, parameter complexity, etc...)

- Note that the meausure the performance (evaluation) with an approach that seems more robust than yours because check in conditions that are more different or something like that.

- I have used this paper to justify some choices we made that are not varied in our paper. Of course, the best option is to change them also, specially model complexity and number of variables, but given that we did not varied them, we select a middle ground, we made choices that seem good in general for most of the Pinus species.
	- the parameter complexity of GLM and GAM. 
		- We used an intermediate level of complexity according to this paper, and this gave the highest levels of accuracy.
		- Note that this relationship was not consistent across all as‐ sessments: TSS of complex parameterizations was highest when considering ‘common’ species (>300 presence observations).
		- But we selected the same parameters for all species, so it seems that intermediate is a good middle ground for the whole genus. It is, indeed, the middle ground, not low, not too high.
		- Of course, the best option would be to vary this parameter, but we did not consider it. And they recomend to be varied in each dataset.
	- the statement according to which avoding a high number of highly correlated predictors is better for SDMs. ALSO the cite indicating the multicolinearity can impact SDM projections and we limit correlation to belo 0.4
		- they find that the accuracy (AUC and TSS) peaks at 10 predictors, with more reducing the accuracy.
			- although performance declines at high numbers of variables were modest
		- multicolinearity of the predictors is not associated with accuracy, but it is with range change. It increases the disagreement among range change between projections, so it should be below 0.7 according to this study in order to avoid its influence on the projections under climate change. It also increase range loss (figure 4).
			- below 0.4 there is not influence on projections disagreement for range change. Figure S6.14.
	- The predictor/presence ratio used in our study. We at least ensured 1/10
		- They also show that a predictor/presence ratio is better for AUC and TSS.
		- although for species with more presences (e.g., 5/500=0.01) the accuracy is worse.
		- As explained below, we tried to use a similar set of predictors for all species, we cannot just use a lot of predictors for the big ones, because then the ratio predictor/presence woudl be much lower than 0.1 for many species. So we should have to reduce the number of predictors for these species, being a great differences between the modeling of big species respect to the rest.
			- in any case, the accuracy in general was very except for very large and very small range species, so it seems we are in middle ground.
		- Also note that the optimal number of variables is dataset specific and they have presence-absence data, not only-presence data.
			- For sam‐ pling designs similar to the one we worked with, one predictor per 10 presence observations may be ideal (Harrell et al., 1998). However, for less well‐designed survey data such as presence‐ only data, finer grains or steeper environmental gradients, more presence observations per predictor may be necessary
			- we have both, well sampled and big distributions and smaller... so we again select a good compromise, specially given the performance we found.
			- and in general, we have presence-only, so maybe our peak is below 10. They say that maybe for only-presence data it is needed more than 10 preseneces per predictor. So it makes sense to at least reach 10, and go up from there as we work with only-presence data.
		- In addition, note that the increase in the number of variables increased the variability across projections for range change, i.e., more disagreement between projections. This can be see in figure S6.14 even for 4-5 variables and above.
	- to support the utility to use different SDM algorithms, because different types can produce different results
		- they show that the predictions of range change are very influenced by the type of algorirthm.
		- RF tends to give positive predictions of range change, while GLM and GAM give negative.
		- there also differences in the variability between projections of range change between models.

- Summary
	- Based on our results and the considerations discussed above, we formulate three recommendations for including model complexity in ensemble simulations of climate change impact on biodiversity using SDMs
		- SDM algorithms and parameterization complexity: SDM algo‐ rithms as well as parameterization complexity have important consequences on projected distributional change and thus both factors should be varied in ensembles. Appropriate levels of parameterization complexity depend on the dataset at hand, and can be constrained based on model performance. We suggest to run SDM algorithms at least at two levels of parameterization complexity. Under computational constraints, this may go at the cost of using many SDM algorithms.
		- Predictor numbers in ensembles: The number of predictors strongly impacts model performance and can affect disagree‐ ment among range loss projections. Our results suggest that op‐ timal performance may be achieved with around 10 predictors, or one predictor per ten presences, if well‐designed survey data and diverse predictors are available. For studies using presence‐only data and/or exclusively climate predictors, this number may well be lower. The strong dependence of model performance on num‐ ber of variables makes it straight‐forward to optimize this factor for the dataset at hand using block cross‐validation.
		- Multicollinearity: In this study, multicollinearity did not strongly affect the performance of model extrapolations, but it distinctly increased projected range loss and the disagreement among range change projections. We recommend keeping absolute Pearson correlation coefficients below 0.7, a boundary recom‐ mended elsewhere (Dormann et al., 2013), and one above which consequences in projections became clearly visible.

## Introduction

- Efficient mitigation of biodiversity loss from global changes requires a thorough understanding of how species’ ranges are organized in space, and how they will shift in the future.

- Two approaches are commonly employed to establish such understanding: statistical species distribution models (SDMs, Guisan & Zimmermann, 2000) and mechanistic models (e.g. Zurell et al., 2016). Projections of spe‐ cies range shifts using mechanistic models are based on explicitly formulated processes that are presumably relevant to the ecology of the target species, while SDM projections extrapolate relationships identified from statistical structures between occurrences and their environment. In principle, projections from mechanistic models may seem preferable as their careful application may harbour a lower risk that relevant processes are insufficiently captured or corrupted by erroneously identified associations (Merow et al., 2014). However, limited understanding of relevant processes and of the ecology of most species, and/or lack of relevant data to describe it sufficiently well prevent their use in many cases (Guisan & Zimmermann, 2000; Thuiller et al., 2008).

- Despite their limitations, statistical SDMs are therefore likely to remain commonly used to project species re‐ sponses to global change. For this reason, IT IS IMPERATIVE TO COMPRE‐ HEND THE IMPLICATIONS OF THE VARIOUS CONCEPTUAL DECISIONS TAKEN AT THE DIFFERENT STEPS OF THE DEVELOPMENT OF SDM PROJECTIONS.

- Implications of decisions in projection design can be quantified
by comparing the outcomes of alternative setups when projected under climate change (aka projection ensembles). Projection en‐ sembles consist of multiple projections generated by systematically varying the settings at the different steps of their development, such as initial conditions, that is, the presence and (pseudo)absence data used for model training, predictor variables, SDM algorithms, pa‐ rameterization complexity, climate models, or emission scenarios. PROJECTION ENSEMBLES ARE PARTICULARLY USEFUL TO QUANTIFY UNCERTAINTY AND TO OBTAIN CONSENSUS PROJECTIONS,  WHICH ARE ARGUABLY SUPERIOR TO SINGLE MODEL PROJECTIONS (Araújo & New, 2007, but see Dormann et al., 2018). Furthermore, if combined with rigorous model valida‐ tion, projection ensembles can help identifying model designs of relatively high quality.

- Compared to other fields, such as econom‐ ics and climate science, projection ensembles were introduced to species distribution modelling relatively recently (Araújo & New, 2007; Thuiller, 2004), but gained popularity since specialized mod‐ elling platforms became available—such as the R‐package ‘biomod2’ (Thuiller, Lafourcade, Engler, & Araújo, 2009).

- However, not all steps in the development of projection ensembles have received the same level of attention. A literature study of 125 recent papers employing SDM projections revealed that the most frequently varied step was the emission scenario (63% of cases), followed by the climate models used to estimate future climatic conditions (48% of cases) (Figure 1a, for further information see Appendix S1). SDM algorithms and ini‐ tial conditions were also frequently varied (35% and 32% of cases, respectively).
	- this is what we did. Sot he SDM algorithms and initial conditions are not very very frequently varied, but they are varied some time...

- Implications of decisions revolving around model com‐ plexity, on the other hand, were typically not explored, and either left to the defaults of the method applied or taken based on more or less well‐grounded heuristics. Yet, the importance of also varying model complexity in projection ensembles has recently been empha‐ sized by several authors (Boria, Olson, Goodman, & Anderson, 2014; Merow et al., 2014; Werkowska, Márquez, Real, & Acevedo, 2017).

- Most SDM algorithms can be tuned to fit models across a sub‐ stantial range of complexity, from ‘under fit’ models that are not flexible enough to capture the detailed species response to the envi‐ ronment to ‘over fit’ models that ascribe signal to noise, which is par‐ ticularly risky when projecting (Merow et al., 2014; Moreno‐Amat et al., 2015). Even when differences in model performance are minor, projections from complex models can strongly differ from those of simple models (Beaumont et al., 2016; Gregr, Palacios, Thompson, & Chan, 2018; Merow et al., 2014).

- However, systematically varying model complexity across different SDM algorithms is not straight‐ forward, as their different setups do not allow for analogous tuning, and universal measures to directly compare complexity are lacking (García‐Callejas & Araújo, 2016).

- We investigate the roles of three aspects related to model complexity: parameterization complexity, number of variables used, and multicollinearity among variables.

- Parameterization complexity involves modifications of a set of pa‐ rameters, adjusting the level of complexity within SDM algorithms. These variations can be based on the flexibility of response curves or the inclusion of interaction terms in regression techniques and tree complexity in tree‐based methods (Merow et al., 2014). ). Varying parameterization complexity has not been employed routinely in the recent literature. Among the 125 papers that we investigated, it was varied only twice (Figure 1a). Instead, algorithms were mostly run with default parameterizations or else with simplifications of the default flexibility (see also Hao, Elith, Guillera‐Arroita, & Lahoz‐ Monfort...).
	- we did NOT do this, because we used algorithsms with different completxity in general, but we did not modify the parameters within each type of algorithm to have a wide range of complexity.

- Model complexity is also affected by the number of predic‐ tor variables considered as well as their multicollinearity. Adding more predictors to a model increases the amount of signal and noise available to SDM algorithms and typically leads to larger numbers of parameters estimated, and thus more complex models (Merow et al., 2014; Werkowska et al., 2017). However, many algorithms include strategies to eliminate parameters that insufficiently improve model fits, which leads to a saturating relationship between number of variables and model complex‐ ity.

- Particularly many parameters may be eliminated for predic‐ tor sets with high levels of multicollinearity, and thus a limited amount of independent information. Multicollinearity may therefore lead to somewhat simpler models. But investigating the effects of multicollinearity is also of interest because it can compromise parameter estimates which are especially problem‐ atic when models are transferred to situations with different multicollinearity regimes (Dormann

- As ecologically important predictors often show significant levels of collinearity, knowing the maximum level of tolerable collinearity is critical. Among the 125 papers we investigated, the median number of variables included was seven, ranging from two to 37 (Figure 1b). 

- Yet, within the same analysis the numbers were typically not var‐ ied (only in 6% of cases), and if they were, then mainly as a con‐ sequence of recombining variable groups (e.g. climate vs. climate and soil variables) and not to study the impact of numbers of variables. Also, multicollinearity levels were only exceptionally varied (2% of cases), and the heuristics used to limit multicol‐ linearity varied greatly (Figure 1c
	- we varied the number of predictors, but as a consequence of the differences in sample size between species. We did not varied the number within the same species to see the impact.

- In this study, we analyzed a comprehensive ensemble of SDM projections and compared uncertainty associated with the com‐ monly varied decision steps in ensembles (SDM algorithm, emis‐ sion scenario, and climate models) with uncertainty originating from parameterization complexity, number of variables, and mul‐ ticollinearity. Furthermore, we investigated the patterns of model performance, projections of distributional change, and disagree‐ ment of projections of distributional change (i.e. variation from replicated predictor sets) along model complexity gradients. Using survey data for 34 tree species across Europe, we fitted and evalu‐ ated more than 100,000 SDMs with two performance metrics, and generated over 800,000 projections of species distribution ranges that we summarized with two metrics of distributional change.

- points addressed:
	- Which are the most important factors affecting the performance of model extrapolations to ‘novel’ (non‐analogous) conditions, and projections of species distributional change?
	- Are the effects of model complexity on model performance and species distributional change in line with the expectations formu‐ lated in Table 1?


## Methods

- Our analyses consisted of three steps. First, we prepared a compre‐ hensive set of environmental variable combinations. We established a pool of 24 climate variables for both, current and future condi‐ tions, and a pool of 16 soil/terrain variables which we assumed to remain constant until 2080. Based on pairwise Pearson correlation coefficients, we defined 100 combinations of numbers of variables and multicollinearity levels, and screened the realm of possible pre‐ dictor sets with roughly equal numbers of climatic and soil/terrain variables for three replicates per combination. Second, we evaluated and projected a large number of SDM fits (Figure 2). For each com‐ bination of predictor set and species, we fitted four SDM algorithms at three levels of parameterization complexity and evaluated their performance. Then, we projected the fitted models to the conditions in 2061–2080 as projected by four climate models for two emission scenarios and assessed projected species distributional changes. Third, we investigated how model complexity affects model perfor‐ mance, projected distributional change, and disagreement between projections of distributional change.

- RF fits were based on 500 trees, and
	- just like we did, they did not modify this parameter.
	- we cannot use this reference to support 500 trees, becuase they do not generate evidence about this.

- We assessed model performance based on two metrics, True Skill Statistic (TSS, Allouche, Tsoar, & Kadmon, 2006), and area under the curve (AUC, Swets, 1988) These metrics were derived from model projections to ‘novel’ conditions using environmental block cross‐validation (Roberts et al., 2017, Appendix S2.1 in Appendix S2). Block cross‐validation is a comparably tough test enforcing projections to conditions that are somewhat more different than our future environmental conditions were from present conditions, both in terms of the covered ranges and correlation structure (Appendix S2.2 in Appendix S2).
	- model performance measure in conditions even more different than those expected under climate change.

- We used analysis of variance (ANOVAs) to assess the relative contributions of the steps in projection development to uncertainty in model performance and projections of distributional change. These analyses were based on the outputs of the primary analyses plus missing value imputations for number of variables × multicollinearity combinations for which no

- Imputation: We summarized model performance and distributional change in the number of variables × multicollinearity space to investigate their patterns and to generate estimates for missing values. For estimat‐ ing missing values, we represented the combinations of number of variables and multicollinearity bins in a 10 × 10 pixel space sepa‐ rately for each combination of species, SDM algorithm, and param‐ eterization complexity (plus emission scenario and climate model in the case of distributional change). Then, we pixel‐wise summarized model performance and projected distributional change estimates from the three replicates by median and approximated pixels with missing values with bilinear interpolations from neighbouring pix‐ els. To investigate patterns, we combined original data and imputed missing values, and similarly summarized pixels by median and inter‐ quartile range (IQR) but for pooled estimates from all species.

- We used ANOVA to quantify the relative contributions of the differ‐ ent sources of uncertainty in projection ensembles. We ran ANOVAs with model performance metrics (TSS and AUC) and with distribu‐ tional change estimates (range loss, range change) as response. For model performance ANOVAs, we compared the contributions of number of variables, multicollinearity, parameterization complex‐ ity, and SDM algorithm. For distributional change ANOVAs, we ad‐ ditionally considered the effects of the two different climate models and the two different emission scenarios. In order to have compara‐ ble level numbers for the different factors, we aggregated number of variables and multicollinearity to two levels: low and high levels of multicollinearity were distinguished by a third quartile of |r| of 0.5 while the group of low numbers of variables included 3–7 and the group of high numbers 8–12. This aggregation resulted in 75 poten‐ tial predictor sets for each of the four combinations of aggregated nvar and multicollinearity levels. To account for non‐independence resulting from the nestedness of SDM algorithm and parameteriza‐ tion complexity, we additionally considered their linear interaction. We accounted for species identity through a random intercept. ANOVAs were based on Bayesian generalized linear mixed models, fitted with the Integrated Nested Laplace Approximations (INLA) ap‐ proach (Rue, Martino, & Chopin, 2009). Instead of p‐values, which are not helpful for large sample sizes, we used parameter uncertainty in the posterior distributions to assess how distinct mean sums of squares of the different factors were. We estimated mean sums of squares 1,000 times based on resampled parameter estimates from the posterior distributions of the fitted INLA models, and report me‐ dians and 95% confidence intervals. For response variables bounded by zero and one (AUC, range loss) we assumed errors to follow a beta distribution, otherwise normal error distribution was assumed
	- I guess they used parameter uncertainty in the posterior distributions to assess how distinct mean sums of squares of the different factors were


## Results

### Model performance

#### Analyses of variance

- I do not fully understand how they calculated the importance of each variable. Because of this, I am not adding notes about this.

#### Analysis of patterns

- Overall, TSS measured under environmentally extrapolating block cross‐validation was highest for parameterizations of intermedi‐ ate complexity, showed a unimodal relationship with number of variables, and no clear relationship with multicollinearity (Figure 4). 

- TSS increase was steep for models built with three to five variables, started levelling‐off for models built on >5 variables, and typically peaked at 10 or 11 variables (Figure 4b). 
	- it is clear that more than 10 predictors is not good. 
	- We used a maximum of 5, but you have to take into account that we used the same sets of variables for many species. to me is rasonable to leave the number of predictors not so high, eve if we do not get the highest highest accuracy.  we tried to use a similar set of predictors for all species, we cannot just use a lot of predictors for the big ones, because then the ratio predictor/presence woudl be much lower than 0.1 for many species. So we should have to reduce the number of predictors for these species, being a great differences between the modeling of big species respect to the rest.
	- Note that
		- the accuracies of our models are really high.
		- only species with very large and very small ranges have low accuracy, so it seems we are in middle ground.

- Fits of intermediate and high complexity achieved notably higher TSS than those of low complex‐ ity. Their TSS was similar if no more than five variables were included; otherwise fits of intermediate parameterization complexity outper‐ formed complex fits (Figure 4b)
	- We have used the intermediate levels of GLM (second order polynomials) and GAM (smooth term of 4, 3 is intermediate level in this paper). In GLM, I think glm from R automatically run the linear component if the poly(2) is used, just like they did in this paper. In RF I am not sure if we changed the number of observations in the last node.
	- The figure 4b and c shows that models with intermediate levels of complexity give the highest TSS. The same is shown for AUC in figure 5.2.3 b and c. In the panels a) of both figures, you can see how the middle column (intermediate) is specially higher respect to simple and complex in the case of GLM and GAM.
	- so it seems we have used a good level of parameter complexity in GLM and GAM.

- Interquartile range of TSS tended to decrease with parameterization complexity, in particular for GLMs and GAMs (Figure S5.5 in Appendix S5), indicating that under these conditions the type of predictor variable used had a comparably low impact on performance.
	- I understand that interquartile range is a measure of the variability across replicates.
	- As the complexity increases, there is less variation between replicates.

- TSS of fits of ‘common’ species was on av‐ erage slightly lower than overall TSS, but the patterns were gener‐ ally similar (Figure S5.6 in Appendix S5). However, ‘common’ species fits showed a weak negative relationship with multicollinearity, and among those fits complex parameterizations achieved highest TSS.
	- Ok, so for common species more complex could be the best, but in general, we can say that intermediate is better.

- Model performance patterns remained similar when assessed by AUC, and when relative instead of absolute predictor numbers were considered. 
	- As TSS, AUC showed a unimodal relationship with num‐ ber of variables with highest scores at 10 or 11 variables, it peaked for parameterizations of intermediate complexity, and showed no clear relationship with multicollinearity.
		- again 10 predictors 
	- Relationships with model performance were also mostly unimodal when the number of variables per presence observation rather than the absolute number of variables was considered (Figure S5.8 in Appendix S5): TSS and AUC were typically highest if one predictor was used for about 10 presence observations.
		- in the figure S5.2.4 you can see that the best AUCs are obtained close to the end of the X axis, which is 0.1 (1 variable/10 predictors)
		- the pattern is very clear to models with intermediate complexity of parameters (our case in GLM and GAM). As we get closer to 1/10=0.1, the TSS and AUC is higher. But for species with more presences, for example 5/500=0.01 we are at lower accuracy.
		- Again, we used tried to use a similar set of predictors for all species, we cannot just use a lot of predictors for the big ones, because then the ratio predictor/presence woudl be much lower than 0.1 for many species. So we should have to reduce the number of predictors for these species, being a great differences between the modeling of big species respect to the rest.
		- only species with very large and very small ranges have low accuracy, so it seems we are in middle ground.
		- at least, we can use this evidence to justify that in small-size species, we ensured at least 10 presences per predictor.

### Species distribution change

#### Analyses of variance

- I do not fully understand how they calculated the importance of each variable. Because of this, I am not adding notes about this.

#### Analysis of patterns

- Higher fractions of ranges were projected to be lost by fits with more complex parameterizations and predictor sets with elevated levels of multicollinearity (Figure 6). 

- On average pa‐ rameterizations of intermediate complexity projected a median range loss that was 16% higher than that of parameterizations of low complexity; fits with parameterizations of high complexity projected another 5% increase (Figure 6b,c). These differences were driven by projections of GLMs and GAMs which were par‐ ticularly affected by parameterization complexity (Figure S6.11 in Appendix S6). Median projected range loss also increased by 10% for predictor sets with a third quartile of r larger than 0.5 (Figure 6b).
	- you can see in figure 6a, more green in the upoper right pannel of GLM and GAM. Also in panels b y d, whith higher median range loss for complex models.
	- multicolinearity was not important for AUC and TSS, but it is important for range loss.
	- therefore, multicolineraty is relevant also for future projections.

- As for range loss, patterns of range change responded to parameterization complexity, with models fitted with parameterizations of low complexity project‐ ing net range gains, and those fitted with other complexity levels projecting net range losses (Figures S6.14 and S6.15 in Appendix S6). However,

- However, SDM algorithms caused even larger differences in projected range change: RF projections estimated notable net range gains (33% on average), whereas all other algorithms over‐ all projected little change or net range losses (3%, −5%, and −8% on average for GLM, GAM and GBM, respectively).
	- you can see this in Figure S6.14 and S6.15. Random forest shows clearly positive values of median range change, while the other models tend to show negative values.
	- Another piece of evidence in support of using different algorithms.

- Responses of projected distributional change were much more pronounced when relative rather than absolute numbers of variables were considered (Figure S6.16 in Appendix S6). Range loss showed a concave relationship with number of variables per presence observation with minima at between 10 and 25 presences per variable (1/25=0.04 and 1/10=0.1); range change projections, in contrast, peaked at these numbers.
	- meaning that there are more clear patterns when using the number of predictors relative to the number of observations.
	- i an not sure if this is preferable or not...

- The interquartile range of range loss projections varied consid‐cerably for different combinations of SDM algorithm and parame‐ terization complexity, and showed a weakly positive relationship with both number of variables and multicollinearity (Figure 6). The IQR of range loss projections was higher for GLMs and GAMs than for GBMs and RF (Figure 6 and Figure S6.11 in Appendix S6). Relationships with multicollinearity and number of variables were both increasing but rather weak, while range loss IQR was slightly higher for parameterizations of intermediate complexity than for those of high or low complexity. Patterns
	- again, the type of the model influence, in this case influences the variability of range loss across projections
	- also there is a weak infleunce by the number of varaibles and multicolinearty. Increasing both the variation across projections.

- Patterns were similar among ‘common’‐species models, although among them simple fits were associated with highest range loss IQR (Figures S6.12 and S6.13 in Appendix S6).

- IQR of projected range change also tended to increase with number of variables and multicollinearity, and decreased with parameterization complexity (Figures S6.14 and S6.15 in Appendix S6). Range change IQR was furthermore especially high for RF projections.
	- infleunce of the type of model and also number of variables and multicollinearity. In this case, more variables and more multicolinearity increases the variability across projections. Respect models, RF has higher variability.


## Discussion

- Read, but added notes along the previous lines.