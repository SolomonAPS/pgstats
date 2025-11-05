# Population Genetics Statistics: Mathematical Formulas

This document provides the complete mathematical formulas for all statistics implemented in pgstats, with proper notation and textbook references.

## Table of Contents

1. [Notation and Definitions](#notation-and-definitions)
2. [Site Frequency Spectrum](#site-frequency-spectrum)
3. [Theta Estimators (Diversity Measures)](#theta-estimators)
4. [Neutrality Tests](#neutrality-tests)
5. [Linkage Disequilibrium Statistics](#linkage-disequilibrium-statistics)
6. [Haplotype Statistics](#haplotype-statistics)
7. [References](#references)

---

## Notation and Definitions

### Basic Quantities

| Symbol | Definition |
|--------|------------|
| $n$ | Sample size (number of chromosomes) |
| $S$ | Number of segregating sites |
| $L$ | Sequence length (callable bases) |
| $\xi_i$ | Number of sites with $i$ derived alleles (unfolded SFS) |
| $\eta_i$ | Number of sites with minor allele count $i$ (folded SFS) |
| $\zeta_1$ | Number of singletons (sites with exactly 1 derived allele) |

### Harmonic Numbers and Coefficients

$$a_1 = \sum_{i=1}^{n-1} \frac{1}{i}$$

$$a_2 = \sum_{i=1}^{n-1} \frac{1}{i^2}$$

$$b_1 = \frac{n+1}{3(n-1)}$$

$$b_2 = \frac{2(n^2 + n + 3)}{9n(n-1)}$$

$$b_n = \sum_{i=1}^{n-1} \frac{1}{i^2} = a_2$$

**Note**: In Walsh & Lynch (2018), $b_n$ refers to the harmonic sum of squares, which is equivalent to $a_2$ in our notation.

### Binomial Coefficient

$$\binom{n}{2} = \frac{n(n-1)}{2}$$

---

## Site Frequency Spectrum

### Unfolded SFS

The unfolded site frequency spectrum requires knowledge of the ancestral state:

$$\xi_i = \text{number of sites with } i \text{ derived alleles, for } i = 1, 2, \ldots, n-1$$

**Reference**: Wakeley (2009), Section 4.3.1

### Folded SFS

When ancestral state is unknown, we use the minor allele frequency:

$$\eta_i = \text{number of sites with minor allele count } i, \text{ for } i = 1, 2, \ldots, \lfloor n/2 \rfloor$$

**Reference**: Wakeley (2009), Section 4.3.1

### Number of Segregating Sites

$$S = \sum_{i=1}^{n-1} \xi_i$$

**Reference**: Wakeley (2009), Equation 4.38

---

## Theta Estimators

Theta ($\theta = 4N_e\mu$) is the population mutation parameter, where $N_e$ is effective population size and $\mu$ is mutation rate per site per generation.

### 1. Theta Pi ($\theta_\pi$) - Nucleotide Diversity

**Based on**: Average pairwise differences

$$\theta_\pi = \frac{1}{\binom{n}{2}} \sum_{i=1}^{n-1} i(n-i)\xi_i$$

**Normalized by sequence length**:

$$\hat{\theta}_\pi = \frac{\theta_\pi}{L}$$

**Reference**: Walsh & Lynch (2018), Equation 9.4; Wakeley (2009), Equation 4.39

**Interpretation**: Expected heterozygosity; average number of pairwise differences per site.

---

### 2. Theta W ($\theta_w$) - Watterson's Estimator

**Based on**: Number of segregating sites

$$\theta_w = \frac{S}{a_1}$$

**Normalized by sequence length**:

$$\hat{\theta}_w = \frac{\theta_w}{L} = \frac{S}{a_1 L}$$

**Reference**: Walsh & Lynch (2018), Equation 9.5; Wakeley (2009), Equation 4.40

**Interpretation**: Unbiased estimator of $\theta$ under neutral Wright-Fisher model.

---

### 3. Theta H ($\theta_h$) - Fay and Wu's Estimator

**Based on**: High-frequency derived alleles

$$\theta_h = \frac{1}{\binom{n}{2}} \sum_{i=1}^{n-1} i^2\xi_i$$

**Normalized by sequence length**:

$$\hat{\theta}_h = \frac{\theta_h}{L}$$

**Reference**: Walsh & Lynch (2018), Equation 9.6; Fay & Wu (2000), Equation 2

**Interpretation**: Emphasizes high-frequency variants; sensitive to recent selective sweeps.

---

### 4. Theta L ($\theta_L$) - Zeng's Estimator

**Based on**: Intermediate-frequency alleles

$$\theta_L = \frac{1}{n-1} \sum_{i=1}^{n-1} i \cdot \xi_i$$

**Normalized by sequence length**:

$$\hat{\theta}_L = \frac{\theta_L}{L}$$

**Reference**: Walsh & Lynch (2018), Equation 9.7; Zeng et al. (2006), Equation 9.28a

**Interpretation**: Weights sites by allele frequency; useful for detecting selective sweeps.

---

## Neutrality Tests

These statistics test for deviations from neutral evolution by comparing different theta estimators.

### 1. Tajima's D

**Compares**: Pairwise differences vs. segregating sites

$$D = \frac{\theta_\pi - \theta_w}{\sqrt{\text{Var}(\theta_\pi - \theta_w)}}$$

**Variance**:

$$\text{Var}(\theta_\pi - \theta_w) = c_1 S + c_2 S(S-1)$$

where:

$$c_1 = b_1 - \frac{1}{a_1}$$

$$c_2 = b_2 - \frac{n+2}{a_1 n} + \frac{a_2}{a_1^2}$$

**Reference**: Walsh & Lynch (2018), Equation 9.8; Tajima (1989); Wakeley (2009), Equation 4.35

**Interpretation**:
- $D = 0$: Consistent with neutral evolution
- $D > 0$: Excess of intermediate-frequency alleles (balancing selection or population contraction)
- $D < 0$: Excess of rare alleles (purifying selection, population expansion, or selective sweep)

---

### 2. Fu and Li's D* (Folded)

**Compares**: Segregating sites vs. singleton frequency

$$D^* = \frac{\frac{S}{a_1} - \frac{n-1}{n}\eta_1}{\sqrt{\text{Var}(D^*)}}$$

**Variance**:

$$\text{Var}(D^*) = \alpha^* S + \beta^* S(S-1)$$

where:

$$c_n = \frac{n+1}{n} - \frac{1}{a_1}$$

$$d_n = \frac{b_2}{a_1^2} - \frac{2}{n}\left(1 + \frac{1}{a_1} + \frac{a_1}{n}\right) - \frac{1}{n^2}$$

$$\beta^* = \frac{1}{a_1^2 + b_2} \cdot d_n$$

$$\alpha^* = \frac{1}{a_1} \cdot c_n - \beta^*$$

**Reference**: Walsh & Lynch (2018), Equation 9.26b; Fu & Li (1993)

**Interpretation**: Tests for an excess or deficit of singleton mutations relative to total segregating sites.

---

### 3. Fu and Li's D (Unfolded)

**Compares**: Segregating sites vs. derived singletons (requires ancestral state)

$$D = \frac{\frac{S}{a_1} - \zeta_1}{\sqrt{\text{Var}(D)}}$$

Uses same variance formula as $D^*$ but with derived singleton count $\zeta_1$.

**Reference**: Walsh & Lynch (2018), Equation 9.26c; Fu & Li (1993)

---

### 4. Fu and Li's F* (Folded)

**Compares**: Pairwise differences vs. singleton frequency

$$F^* = \frac{\theta_\pi - \frac{n-1}{n}\eta_1}{\sqrt{\text{Var}(F^*)}}$$

**Variance**:

$$\text{Var}(F^*) = \alpha_F S + \beta_F S(S-1)$$

where:

$$\beta_F = \frac{1}{a_1^2 + b_2} \left[\frac{2n^3 + 110n^2 - 255n + 153}{9n^2(n-1)} + \frac{2(n-1)a_1}{n^2} - \frac{8b_2}{n}\right]$$

$$\alpha_F = \frac{1}{a_1} \left[\frac{4n^2 + 19n + 3 - 12(n+1)a_{n+1}}{3n(n-1)}\right] - \beta_F$$

**Reference**: Walsh & Lynch (2018), Equation 9.26e; Fu & Li (1993)

**Interpretation**: Combines information from pairwise differences and singleton frequency.

---

### 5. Fu and Li's F (Unfolded)

**Compares**: Pairwise differences vs. derived singletons

$$F = \frac{\theta_\pi - \zeta_1}{\sqrt{\text{Var}(F)}}$$

Uses same variance formula as $F^*$ but with derived singleton count $\zeta_1$.

**Reference**: Walsh & Lynch (2018), Equation 9.26f; Fu & Li (1993)

---

### 6. Zeng's E

**Compares**: Zeng's theta vs. Watterson's theta

$$E = \frac{\theta_L - \theta_w}{\sqrt{\text{Var}(\theta_L - \theta_w)}}$$

**Variance** (using scaled $\theta$ from Equation 9.21b):

$$\theta = \frac{S}{a_1}, \quad \theta^2 = \frac{S(S-1)}{a_1^2 + b_1}$$

$$\text{Var}(\theta_L - \theta_w) = \left[\frac{n}{2(n-1)} - \frac{1}{a_n}\right]\theta + \left[\frac{b_n}{a_n^2} + 2\left(\frac{n}{n-1}\right)^2 b_n - \frac{2(nb_n - n + 1)}{(n-1)a_n} - \frac{3n+1}{n-1}\right]\theta^2$$

**Reference**: Walsh & Lynch (2018), Equation 9.28c; Zeng et al. (2006)

**Interpretation**: Sensitive to recent selective sweeps; negative values suggest recent positive selection.

---

### 7. Fay and Wu's H

**Compares**: Pairwise differences vs. high-frequency alleles

$$H = \frac{\theta_\pi - \theta_h}{\sqrt{\text{Var}(\theta_\pi - \theta_h)}}$$

**Variance** (using scaled $\theta$ from Equation 9.21b):

$$\theta = \frac{S}{a_1}, \quad \theta^2 = \frac{S(S-1)}{a_1^2 + b_1}$$

$$\text{Var}(\theta_\pi - \theta_h) = u_H \theta + v_H \theta^2$$

where:

$$u_H = \frac{n-2}{6(n-1)}$$

$$v_H = \frac{18n^2(3n+2)b_{n+1} - (88n^3 + 9n^2 - 13n + 6)}{9n(n-1)^2}$$

and $b_{n+1} = \sum_{i=1}^{n} \frac{1}{i^2} = a_2(n+1)$ is the harmonic sum of squares up to $n$.

**Reference**: Walsh & Lynch (2018), Equation 9.27b; Fay & Wu (2000); Zeng et al. (2006)

**Interpretation**: Negative values indicate an excess of high-frequency derived alleles, suggesting recent positive selection.

---

## Linkage Disequilibrium Statistics

### 1. D (Coefficient of Linkage Disequilibrium)

For two loci A and B with alleles $A_1, A_2$ and $B_1, B_2$:

$$D = P(A_1B_1) - P(A_1)P(B_1)$$

where $P(A_1B_1)$ is the frequency of haplotype $A_1B_1$, and $P(A_1), P(B_1)$ are allele frequencies.

**For unphased genotype data**, we use the correlation-based approach:

$$D = r \sqrt{p_A(1-p_A) \cdot p_B(1-p_B)}$$

where $r$ is the correlation coefficient between genotypes at the two loci.

**Reference**: Walsh & Lynch (2018), Equation 9.1

---

### 2. D' (Standardized D)

$$D' = \frac{D}{D_{\max}}$$

where:

$$D_{\max} = \begin{cases}
\min(p_A p_b, p_a p_B) & \text{if } D > 0 \\
\min(p_A p_B, p_a p_b) & \text{if } D < 0
\end{cases}$$

**Range**: $-1 \leq D' \leq 1$

**Reference**: Walsh & Lynch (2018), Equation 9.2

**Interpretation**: $|D'| = 1$ indicates complete LD (no recombination between loci).

---

### 3. r² (Squared Correlation Coefficient)

$$r^2 = \frac{D^2}{p_A(1-p_A) \cdot p_B(1-p_B)}$$

**Range**: $0 \leq r^2 \leq 1$

**Reference**: Walsh & Lynch (2018), Equation 9.3

**Interpretation**: Measures the proportion of variance at one locus explained by the other; $r^2 = 1$ indicates perfect LD.

---

### 4. Omega Statistic ($\omega$)

**Detects**: Recombination breakpoints

For a set of SNPs, partition into left (L) and right (R) regions at each potential breakpoint $k$:

$$\omega_k = \frac{|\bar{r}^2_L - \bar{r}^2_R|}{|\bar{r}^2_{LR}|}$$

where:
- $\bar{r}^2_L$ = mean $r^2$ within left region
- $\bar{r}^2_R$ = mean $r^2$ within right region  
- $\bar{r}^2_{LR}$ = mean $r^2$ between left and right regions

The omega statistic is:

$$\omega = \max_k \omega_k$$

**Reference**: Kim & Nielsen (2004)

**Interpretation**: High $\omega$ values indicate a recombination breakpoint, suggesting a selective sweep.

---

## Haplotype Statistics

### 1. Haplotype Diversity (H)

**Definition**: Probability that two randomly chosen haplotypes differ

$$H = 1 - \sum_{i=1}^{k} p_i^2$$

where $p_i$ is the frequency of the $i$-th haplotype and $k$ is the number of distinct haplotypes.

**Reference**: Walsh & Lynch (2018)

**Interpretation**: Similar to heterozygosity but for haplotypes; ranges from 0 (no diversity) to ~1 (high diversity).

---

### 2. Garud's H1

**Definition**: Homozygosity of the most common haplotype

$$H_1 = \sum_{i=1}^{k} p_i^2$$

where haplotypes are ordered by frequency: $p_1 \geq p_2 \geq \ldots \geq p_k$.

**Reference**: Garud et al. (2015)

**Interpretation**: High $H_1$ indicates low haplotype diversity; sensitive to hard selective sweeps.

---

### 3. Garud's H12

**Definition**: Combined frequency of the two most common haplotypes

$$H_{12} = (p_1 + p_2)^2 + \sum_{i=3}^{k} p_i^2$$

**Reference**: Garud et al. (2015)

**Interpretation**: Distinguishes hard sweeps (high $H_1$) from soft sweeps (high $H_{12}$ but moderate $H_1$).

---

### 4. Garud's H123

**Definition**: Combined frequency of the three most common haplotypes

$$H_{123} = (p_1 + p_2 + p_3)^2 + \sum_{i=4}^{k} p_i^2$$

**Reference**: Garud et al. (2015)

**Interpretation**: Further extends the ability to detect soft sweeps from multiple haplotypes.

---

### 5. Garud's H2/H1

**Definition**: Ratio of homozygosity excluding the most common haplotype

$$\frac{H_2}{H_1} = \frac{\sum_{i=2}^{k} p_i^2}{\sum_{i=1}^{k} p_i^2}$$

**Reference**: Garud et al. (2015)

**Interpretation**: 
- Low $H_2/H_1$: Single dominant haplotype (hard sweep)
- High $H_2/H_1$: Multiple common haplotypes (soft sweep or balancing selection)

---

## References

### Textbooks

1. **Walsh, B. & Lynch, M. (2018)**. *Evolution and Selection of Quantitative Traits*. Oxford University Press.
   - Chapter 9: Equations for theta estimators and neutrality tests

2. **Wakeley, J. (2009)**. *Coalescent Theory: An Introduction*. Roberts & Company Publishers.
   - Chapter 4: Site frequency spectrum and theta estimators

### Papers

3. **Tajima, F. (1989)**. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics*, 123(3), 585-595.
   - Original Tajima's D statistic

4. **Fu, Y. X., & Li, W. H. (1993)**. Statistical tests of neutrality of mutations. *Genetics*, 133(3), 693-709.
   - Fu and Li's D, D*, F, and F* statistics

5. **Fay, J. C., & Wu, C. I. (2000)**. Hitchhiking under positive Darwinian selection. *Genetics*, 155(3), 1405-1413.
   - Fay and Wu's H statistic and theta_h estimator

6. **Kim, Y., & Nielsen, R. (2004)**. Linkage disequilibrium as a signature of selective sweeps. *Genetics*, 167(3), 1513-1524.
   - Omega statistic for detecting recombination breakpoints

7. **Zeng, K., Fu, Y. X., Shi, S., & Wu, C. I. (2006)**. Statistical tests for detecting positive selection by utilizing high-frequency variants. *Genetics*, 174(3), 1431-1439.
   - Zeng's E statistic and theta_L estimator

8. **Garud, N. R., Messer, P. W., Buzbas, E. O., & Petrov, D. A. (2015)**. Recent selective sweeps in North American Drosophila melanogaster show signatures of soft sweeps. *PLoS Genetics*, 11(2), e1005004.
   - Garud's H1, H12, H123, and H2/H1 statistics

---

## Implementation Notes

### Missing Data

All statistics properly handle missing data:

1. **Per-site sample sizes**: Theta estimators use the actual number of non-missing samples at each site
2. **Variance calculations**: Use the maximum sample size across variants for neutrality test variances (following sgkit and scikit-allel conventions)

### Callable Sites

Theta estimators are normalized by the **callable sequence length** $L$:

$$\hat{\theta} = \frac{\theta}{L}$$

where $L$ excludes regions that cannot be genotyped (e.g., low coverage, repetitive sequences). This is specified via BED files in the CLI.

### Folded vs. Unfolded

- **Unfolded statistics** (D, F, H, E) require knowledge of the ancestral state
- **Folded statistics** (D*, F*) use the minor allele frequency and don't require ancestral state information
- By default, Fu and Li's statistics use the folded version (D*, F*)

