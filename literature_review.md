
# A Literature-Based Interrogation of the Microbial Opsin Structural Landscape: Addressing Gaps and Uncovering Novel Functional Frontiers

## Table 1: Summary Status of Missing Function-Domain Combinations

| Molecular Function | Domain | Status | Key Evidence / Rationale |
|-------------------|---------|---------|--------------------------|
| Anion Channel | Archaea | True Absence | Channel function is a well-documented eukaryotic innovation for rapid photosignaling; Archaea and Bacteria utilize a distinct two-component sensory-pump signaling architecture. |
| Cation Channel | Archaea | True Absence | Channel function is a well-documented eukaryotic innovation for rapid photosignaling; Archaea and Bacteria utilize a distinct two-component sensory-pump signaling architecture. |
| Cation Channel | Bacteria | True Absence | Channel function is a well-documented eukaryotic innovation for rapid photosignaling; Archaea and Bacteria utilize a distinct two-component sensory-pump signaling architecture. |
| Chloride Pump | Archaea | Fillable Gap | Archetypal Halorhodopsins (HRs) are well-characterized archaeal Cl⁻ pumps, first discovered in Halobacterium salinarum. |
| Sodium Pump | Archaea | True Absence | Light-driven Na⁺-pumping is a specialized bacterial innovation defined by the NDQ motif; extensive genomic surveys show no evidence of this function in Archaea. |
| Sodium Pump | Eukaryota | True Absence | Light-driven Na⁺-pumping is a specialized bacterial innovation defined by the NDQ motif; extensive genomic surveys show no evidence of this function in Eukaryota. |
| Unknown | Archaea | Fillable Gap | Numerous candidates exist, including Heliorhodopsins, newly discovered Heimdallarchaeial rhodopsins, CryoRhodopsins, and the broad class of Alt-Rhodopsins. |
| Unknown | Eukaryota | Fillable Gap | Numerous candidates exist, including the newly discovered Apusomonad rhodopsins, Xenopsins, orphan fungal and animal opsins, and the broad class of Alt-Rhodopsins. |


## I. Introduction: From Known Structures to Known Unknowns

### 1.1. Contextualizing the Analysis

The provided structural analysis of the microbial opsin family represents a significant step forward in mapping the functional landscape of these remarkable molecular machines. By integrating a foundational set of 62 experimental structures with 57 non-redundant structures derived from high-quality computational prediction, the analysis effectively doubles the available structural data. This approach is emblematic of a paradigm shift in modern structural biology, where computational methods are no longer merely adjuncts to experimental techniques but are indispensable tools for navigating the vast, uncharted territories of protein sequence space.1 The finding that these predictions fill eight critical gaps where no experimental structures existed underscores the transformative power of this synergy. The resulting combined dataset of 119 structures provides a far more comprehensive, albeit still incomplete, view of the opsin superfamily.

### 1.2. Report Objective and Scope

This report provides a thorough literature-based interrogation of the remaining gaps identified in the aforementioned analysis. The primary objective is to systematically evaluate the eight function-domain combinations for which no structural representative—experimental or predicted—currently exists. For each identified gap, this review will determine whether it represents:

1. **A "fillable gap" or "gap in data,"** where robust evidence confirms the biological existence of the function-domain pair, but structural characterization is lacking. For these, specific, well-documented protein examples will be provided as prime candidates for future study.

2. **A "true absence" or "gap in nature,"** where the function-domain combination is not known to exist. For these, the underlying evolutionary, biophysical, and functional rationale for its absence will be elucidated based on the current body of scientific literature.
Furthermore, this report will venture into the frontiers of opsin discovery by systematically exploring the "Unknown" functional category. The explosion of genomic and metagenomic sequencing has led to the identification of numerous novel opsin families whose functions are either partially or completely uncharacterized.3 These proteins, often referred to as "orphan" opsins, represent the most promising candidates for expanding our understanding of opsin functional diversity. This review will identify, describe, and contextualize these novel families, providing a roadmap for the next wave of structural and functional characterization.

## II. Analysis of Gaps in Ion Pumping Functions

Light-driven ion pumps are the canonical function of microbial rhodopsins, serving as single-component engines for converting light energy into electrochemical potential.7 The analysis of gaps within this functional class reveals both simple data omissions and profound evolutionary patterns.

### 2.1. The Archaeal Chloride Pump: Re-integrating a Foundational Discovery

The initial analysis identifies a missing entry for the (Chloride Pump x Archaea) combination. A review of the literature demonstrates unequivocally that this is a fillable gap. The light-driven, inward-pumping chloride pump is not only known in Archaea but is one of the foundational discoveries in the field, second only to bacteriorhodopsin.
The archetypal archaeal chloride pump is **Halorhodopsin (HR)**, discovered in the extremely halophilic archaeon *Halobacterium salinarum* in 1977.⁸ Its function as a light-driven, inward-directed electrogenic chloride (Cl⁻) pump is firmly established through decades of biophysical and structural studies.⁷ The physiological role of HR is primarily to maintain osmotic balance by accumulating intracellular Cl⁻ to counteract the high external salt concentrations of the hypersaline environments inhabited by haloarchaea.⁸
While *H. salinarum* HR (HsHR) is the original and most-studied example, another key candidate for characterization is Halorhodopsin from *Natronomonas pharaonis* (NpHR). NpHR is a well-characterized homologue that has become a workhorse in the field of optogenetics, where it is widely used as a tool for robust, light-induced neural silencing.¹⁰ Its advantageous properties, including facile heterologous expression and purification, make it an excellent target for structural studies.
The absence of this function-domain pair from a limited dataset can be explained by the "patchy" distribution of rhodopsin genes among haloarchaeal lineages.¹³ This distribution pattern is not indicative of rarity but rather reflects a complex evolutionary history dominated by both Lateral Gene Transfer (LGT) and frequent instances of gene loss.¹³ Phylogenetic analyses suggest that the four major rhodopsin functions (proton pump, chloride pump, and two sensory types) were likely present in the last common ancestor of haloarchaea, with subsequent, repeated losses in many descendant lineages.¹³
The existence of archaeal HRs also provides a powerful example of convergent evolution when compared with their bacterial counterparts. Archaeal HRs are defined by a key three-residue sequence in their active site, the **TSA motif** (Threonine-Serine-Alanine).¹⁵ In 2014, a functionally analogous but evolutionarily distinct class of chloride pumps was discovered in marine flavobacteria, exemplified by *Nonlabens marinus* rhodopsin-3 (NM-R3).¹² Phylogenetic analysis confirms that these bacterial pumps belong to a "distinct phylogenetic lineage quite distant from archaeal inward Cl⁻-pumping rhodopsins".¹² This functional convergence is reflected at the molecular level; the bacterial pumps achieve chloride transport using a different active site solution, the **NTQ motif** (Asparagine-Threonine-Glutamine).¹² More recently, a third group of cyanobacterial chloride pumps has been identified, represented by *Mastigocladopsis repens* rhodopsin (MastR), which utilizes yet another sequence signature, the **TSD motif** (Threonine-Serine-Aspartate).¹¹ The independent evolution of the same biochemical function at least three separate times, each time employing a different set of critical residues, highlights the remarkable functional plasticity of the seven-transmembrane (7TM) rhodopsin scaffold.
**In conclusion, the absence of an archaeal chloride pump in the analyzed dataset is a sampling artifact. This is a fillable data gap, and HsHR or NpHR are the prime, well-documented candidates to include.**

### 2.2. Sodium-Pumping Rhodopsins: A Specialized Bacterial Innovation

The analysis identifies two missing entries related to sodium pumps: (Sodium Pump x Archaea) and (Sodium Pump x Eukaryota). In contrast to the archaeal chloride pump, the literature provides compelling evidence that these are true absences. Light-driven sodium pumping appears to be a specialized function that, based on all current knowledge, has evolved and is restricted to the domain Bacteria.
The first light-driven sodium (Na⁺) pump (NaR), *Krokinobacter eikastus* rhodopsin 2 (KR2), was discovered in 2013 in the marine flavobacterium *K. eikastus*.¹⁵ This was a landmark finding because it overturned the long-held assumption that the positively charged retinal Schiff base at the core of the protein would electrostatically repel the transport of any cation other than a proton (H⁺).¹⁵ KR2 and its homologs function as outward-directed Na⁺ pumps, converting light energy into a sodium motive force that can power other cellular processes, representing a novel form of photoheterotrophy.¹²
The molecular key to this unique function is a distinctive amino acid signature in the third transmembrane helix: the **NDQ motif** (Asparagine-Aspartate-Glutamine).¹² This motif unequivocally distinguishes NaRs from all other known ion-pumping rhodopsins, including proton pumps (which typically have DTD or DTE motifs) and the various families of chloride pumps (TSA, NTQ, or TSD motifs).¹⁵
An exhaustive review of the provided literature reveals no evidence for native, light-driven sodium-pumping rhodopsins in either Archaea or Eukaryota.
**Absence in Archaea:** Discussions of archaeal rhodopsin pumps are consistently focused on the outward proton pump Bacteriorhodopsin (BR) and the inward chloride pump Halorhodopsin (HR).⁷ Even the most recent metagenomic surveys that have uncovered novel archaeal rhodopsin families, such as those from Asgard archaea (Heimdallarchaeial rhodopsins) or from cryo-environments (CryoRhodopsins), have identified only proton pumps or proteins with sensory-like functions.⁶ The NDQ motif characteristic of sodium pumps has not been reported in any archaeal lineage.
**Absence in Eukaryota:** Known eukaryotic rhodopsin pumps are almost exclusively outward proton pumps, identified in fungi (e.g., *Leptosphaeria maculans* rhodopsin, LR) and various phytoplankton such as diatoms.² No naturally occurring eukaryotic sodium pump has been described.
The phylogenetic constraint of NaRs to a specific bacterial clade suggests it is either a relatively recent evolutionary innovation that has not yet disseminated via LGT, or it is a specialized adaptation to a niche whose requirements are not broadly applicable to archaea or eukaryotes. While LGT is a major force in opsin evolution¹³, the specific and complex set of mutations required to re-engineer a proton pump into a sodium pump—solving the fundamental problem of transporting a cation past the positively charged retinal Schiff base—may represent a high evolutionary barrier. This barrier appears to have been surmounted only once, and the resulting functional module has so far remained confined to its bacterial lineage of origin.
**Therefore, the absence of sodium-pumping rhodopsins in Archaea and Eukaryota within the dataset reflects a genuine biological reality. These are true gaps in nature, and this function should be considered a hallmark of bacterial rhodopsin evolution.**

## III. The Great Divide: Evolutionary Divergence of Ion Channels and Pumps

The most significant pattern of absence in the structural dataset is the complete lack of native light-gated channels in the prokaryotic domains (Archaea and Bacteria). This is not an artifact of sampling but rather reflects a fundamental fork in the evolutionary road, where prokaryotes and eukaryotes adopted profoundly different strategies for microbial light-sensing and response.

### 3.1. Biophysical and Evolutionary Principles: The "One-Gate vs. Two-Gates" Model

To comprehend the starkly segregated distribution of channels and pumps across the domains of life, one must first appreciate their fundamental mechanistic and functional opposition. The literature provides a powerful conceptual framework for this distinction: the **"one-gate-versus-two-gates" model** of membrane transport.³²
**A light-driven pump** is an active transporter that operates via an alternating access mechanism.¹⁵ It can be conceptualized as a pore with at least two gates—one facing the extracellular side and one facing the cytoplasmic side—that are strictly coordinated to never be open at the same time. The transport cycle involves the binding of an ion on one side, closure of the first gate, an energy-dependent conformational change in the protein (driven by light absorption), opening of the second gate, and subsequent release of the ion. This process is necessarily slow and deliberate, with transport rates limited by the speed of the protein's conformational changes (on the order of hundreds of events per second).³² This mechanism allows for thermodynamically uphill transport, building and maintaining an ion gradient. A critical feature of this mechanism is the existence of "occluded states," where both gates are closed, trapping the ion within the protein. This fail-safe is essential to prevent the channel-like leakage of ions that would dissipate the very gradient the pump works to build.³²
In stark contrast, **a passive ion channel** is a conduit with a single effective gate. When this gate opens, it forms a continuous, water-filled pathway across the membrane. This allows for extremely rapid, thermodynamically downhill diffusion of ions along their pre-existing electrochemical gradient.³² Ion flux through a single open channel can reach rates of millions to tens of millions of ions per second, orders of magnitude faster than a pump.³² The function of a channel is not to build potential energy in the form of an ion gradient, but to rapidly dissipate that energy for the purpose of signaling, such as triggering a neuronal action potential or a cellular phototactic response.³⁴
The biophysical roles of pumps and channels are thus diametrically opposed. Pumps are metabolic engines that slowly and methodically build electrochemical potential. Channels are rapid signaling switches that consume that potential for fast communication. A pump that leaks is functionally useless; a channel that is slow is equally so. This profound functional dichotomy means that pumps and channels evolved under entirely different selective pressures. The evolutionary pathways to optimize for these two opposing functions are necessarily divergent, and the transition from one to the other would require a significant re-engineering of the protein's core gating mechanism. This provides the theoretical foundation for the clean functional split observed across the domains of life.

### 3.2. Channelrhodopsins: A Eukaryotic Solution for Photosensory Signaling

The structural analysis correctly identifies gaps for (Cation Channel × Archaea), (Cation Channel × Bacteria), and (Anion Channel × Archaea). The scientific literature is unequivocal in confirming that these are true absences. Naturally occurring light-gated ion channels, collectively known as channelrhodopsins (ChRs), are a eukaryotic innovation, found to date exclusively in various lineages of algae and other protists.³⁴
The discovery and function of ChRs are well-documented:
**Cation Channelrhodopsins (CCRs):** The first channelrhodopsins, ChR1 and ChR2, were identified in the green alga *Chlamydomonas reinhardtii*.³⁴ They function as light-gated, non-selective cation channels. Upon illumination with blue light, they open to allow a rapid influx of cations—including H⁺, Na⁺, K⁺, and Ca²⁺—which depolarizes the cell membrane.³⁵ This depolarization serves as the primary electrical signal that guides the organism's phototactic movements.³⁵
**Anion Channelrhodopsins (ACRs):** A second major family of ChRs that are selective for anions was later discovered in cryptophyte algae, with *Guillardia theta* ACR1 (GtACR1) as the prototype.⁴⁵ These channels mediate a light-gated influx of Cl⁻, which hyperpolarizes the membrane of most animal cells. This property has made them exceptionally powerful and widely used tools for the light-based inhibition of neuronal activity in optogenetics.³⁴
The evolution of channel function is itself a story of convergence within the eukaryotic domain. For example, some cation channels found in cryptophytes, termed "bacteriorhodopsin-like cation channelrhodopsins" (BCCRs), show higher sequence homology to archaeal proton pumps than to the canonical CCRs from green algae, suggesting they evolved their channel function independently.³⁶ This pattern of convergent evolution has been dramatically reinforced by the very recent (2025) discovery of **Apusomonad rhodopsins (ApuRs)**, an entirely new family of anion-selective channels from apusomonad protists that are phylogenetically distant from all previously known ChRs.⁵⁰
A critical point of clarification is the term "bacterial channelrhodopsin," which occasionally appears in the literature. A careful reading reveals that this phrase almost invariably refers to an algal channelrhodopsin that has been heterologously expressed in bacteria (such as *E. coli*) for the purpose of protein production and purification for structural or biophysical studies, or in the context of optogenetic applications where the term "bacterial" is used loosely to mean "microbial".⁵² Foundational studies consistently trace the natural origin of all known ChRs to eukaryotic microorganisms.³⁴ This common misnomer is a potential source of confusion but does not constitute evidence for native channel function in prokaryotes.

### 3.3. Prokaryotic Photosensing: The Sensory-Pump Alternative

Prokaryotes are not devoid of sophisticated photosensory capabilities; they simply achieve them through a different evolutionary strategy that leverages the ancestral pump architecture. The prokaryotic solution involves Type 1 rhodopsins that are structurally pump-like but have been repurposed to function as sensors. These are the **Sensory Rhodopsins (SRI and SRII)**, first characterized in *H. salinarum*.⁷
Instead of forming an ion-conducting pore, these sensory rhodopsins function as the photoreceptive subunit of a two-component signaling system. Upon light activation, the SR undergoes a conformational change that allows it to interact with a cognate membrane-embedded transducer protein (Htr).⁸ This protein-protein interaction initiates a cytoplasmic signal transduction cascade, which typically culminates in the modulation of flagellar motor rotation to direct phototactic behavior.¹ This mechanism represents a more conservative evolutionary path than evolving a channel *de novo*. It co-opts an existing, stable 7TM protein fold for a new signaling purpose, essentially functioning as a "broken" pump where the ion translocation pathway is incomplete but the light-driven conformational changes are preserved and repurposed for signal transmission.
**In conclusion, the absence of light-gated ion channels in Archaea and Bacteria is a true biological absence that reflects a deep evolutionary divergence in photosensory signaling strategies. Prokaryotes utilize two-component sensory-pump systems for phototaxis, while eukaryotes evolved single-component, high-throughput ion channels to achieve the same goal. The structural dataset correctly mirrors this fundamental biological reality.**

## IV. Charting the "Unknown": Novel and Uncharacterized Opsin Families

The "Unknown" function category in the structural analysis represents the most dynamic and exciting frontier in opsin research. This is where the vast, unexplored sequence space being uncovered by genomic and metagenomic surveys is yielding novel protein families with functions that challenge our classical understanding of rhodopsins. These "orphan" opsins are the prime candidates for filling the functional and structural gaps in the opsin landscape.
### Table 2: Candidate Opsin Families for the 'Unknown' Functional Category

| Family Name | Domain(s) | Proposed/Known Function | Key Structural/Motif Features | Key References |
|------------|-----------|-------------------------|------------------------------|----------------|
| Heliorhodopsins (HeRs) | Bacteria, Archaea, Eukarya, Viruses | Sensory/Regulatory; Enzyme regulation (GS, photolyase); H⁺ transport (viral only) | Inverted membrane topology (N-in, C-out); ESL motif; Extremely slow photocycle | 54 |
| Alt-Rhodopsins (AltRs) | Bacteria, Archaea, Eukarya, Viruses | Highly diverse; includes K⁺ pumps, Anion Channels, many with unknown function | Substitution at highly conserved Arg82 (BR numbering) by H, K, Q, A, P, S, Y, E, M, W, etc. | 3 |
| Heimdallarchaeial Rs | Archaea (Asgard) | H⁺ pump with light-harvesting carotenoid antenna | Fenestrated structure (Gly at pos. 156); Binds hydroxylated carotenoids (e.g., fucoxanthin) | 28 |
| CryoRhodopsins (CryoRs) | Archaea, Bacteria | Dual-function: inward H⁺ pump / photosensor; UV-switchable | Unique buried Arginine residue; Extremely slow, bimodal photocycle | 6 |
| Apusomonad Rs (ApuRs) | Eukarya (Protists) | Anion-selective channel | Evolved convergently from other ChRs; UV/violet-shifted absorption spectra | 50 |
| Xenopsins | Eukarya (Protostomes) | Visual/Non-visual photoreception; function is poorly characterized | Phylogenetically distinct opsin subgroup; can be co-expressed with r-opsins in mollusks | 59 |
| Orphan Fungal Opsins | Eukarya (Fungi) | Mostly unknown; some implicated in sexual reproduction (e.g., Neurospora Opsin-1) | Large, diverse group with many uncharacterized members beyond known proton pumps | 9 |


### 4.1. The Heliorhodopsin (HeR) Enigma: A New Family with Novel Functions

Heliorhodopsins (HeRs) represent one of the most significant recent discoveries in the field. First identified through functional metagenomics in 2018, they form a phylogenetically distinct family found across all domains of life and even in giant viruses.⁵⁴ Their most revolutionary feature is an **inverted membrane topology**, with the N-terminus facing the cytoplasm and the C-terminus facing the extracellular space—the exact opposite of all previously known Type 1 and Type 2 rhodopsins.⁶³ Initial characterizations of prokaryotic HeRs revealed no ion-pumping activity and extremely slow photocycles, leading to the hypothesis that they function as photosensors.⁵⁵
However, subsequent research has revealed that HeRs are not merely sensors but a new class of light-dependent regulatory proteins, moving them from the "Unknown" category into a new functional space.

- An HeR from *Actinobacteria bacterium* IMCC26103 was shown to physically bind to and regulate the activity of glutamine synthetase, a key enzyme in nitrogen metabolism.⁵⁷
- An HeR from the fungus *Thermomyces flocculiformis* was demonstrated to bind to photolyase, a DNA repair enzyme, and enhance its activity in response to green light, broadening its effective spectrum.⁵⁶
- Intriguingly, a viral HeR from *Emiliania huxleyi* virus 202 was recently shown to function as a light-activated proton transporter, a function not yet observed in any cellular HeR.⁵⁴
HeRs are therefore prime candidates for populating the "Unknown" functional space in Bacteria and Archaea. Their characterization has unveiled an entirely new functional paradigm for microbial opsins that extends beyond simple ion transport or signal transduction via Htr proteins. They can act as light-switchable allosteric modulators of other crucial cellular enzymes, linking environmental light cues directly to metabolic and repair pathways.

### 4.2. Expanding the Blueprint: Alt-Rhodopsins (AltRs) and the "Dark Matter" of the Opsin Universe

The term "Alt-Rhodopsins" (AltRs) was recently proposed as an umbrella classification for any rhodopsin that features a substitution at the highly conserved and functionally critical Arginine 82 residue (using bacteriorhodopsin numbering).³ In canonical proton pumps like bacteriorhodopsin, Arg82 is a key component of the "proton release group," facilitating the transfer of a proton to the extracellular space.⁴
For decades, this residue was considered nearly immutable. However, large-scale genomic and metagenomic surveys have revealed that substitutions at this position are far more common than previously imagined. Residues including Histidine, Lysine, Glutamine, Alanine, Proline, Serine, Tyrosine, Glutamate, and Methionine have all been identified at this position in naturally occurring rhodopsins.³ These AltRs are not confined to a single lineage but are found scattered across all domains of life and viruses.⁴
The functional implications of this variation are profound and place AltRs squarely in the "Unknown" category. While the function of most AltRs remains obscure, the few that have been characterized represent some of the most exciting recent discoveries in the field. For example, the inward-pumping potassium rhodopsins (KCRs) possess a Tryptophan (W) at this position, and some anion channelrhodopsins (ACRs) have a Lysine (K).³ This demonstrates that variation at this single, critical position is a major hotspot for evolutionary innovation, enabling the rhodopsin scaffold to acquire entirely new functions. AltRs are not a single family but rather a principle of variation, representing the vast, uncharacterized "dark matter" of the opsin universe. A clear strategy for discovering novel opsin functions is to systematically mine genomic data for AltR sequences, prioritizing those with unique substitutions from under-sampled phylogenetic clades for functional and structural characterization.

### 4.3. New Frontiers in Archaea: Filling the "Unknown" Gap

Recent metagenomic discoveries have unveiled novel archaeal opsin families that are excellent candidates for the "Unknown" category, each with unique and fascinating properties.
**Heimdallarchaeial Rhodopsins (HeimdallRs):** Discovered in metagenomic assemblies from Asgard archaea—the closest known prokaryotic relatives of eukaryotes—these proteins are rewriting our understanding of the metabolic capabilities of our distant ancestors.²⁸ Functionally, they are proton pumps, but they possess a novel feature: a structural "fenestration" or lateral opening that allows them to non-covalently bind hydroxylated carotenoids (xanthophylls) such as fucoxanthin and lutein.²⁸ These carotenoids act as light-harvesting antennas, absorbing light at wavelengths where retinal does not and transferring the energy to the retinal chromophore, thereby broadening the spectrum of light that can be used for energy generation. This antenna function was previously known only in a few specialized bacterial rhodopsins (xanthorhodopsins). The discovery of this capability in Asgard archaea has profound implications for the evolution of light-harvesting and suggests that the archaeal host of the proto-mitochondrial endosymbiont may have been photoheterotrophic.²
**CryoRhodopsins (CryoRs):** A novel group of rhodopsins recently identified from cold environments like glaciers and cold-adapted organisms.⁶ The characterized CryoR1 exhibits an unprecedented dual functionality, capable of switching between inward proton translocation (a pump function) and a photosensory mode. This functional bimodality is modulated by UV light and is linked to an extremely slow photocycle, which is stabilized by a unique, buried arginine residue not seen in other rhodopsins.
These new archaeal families are not minor variations on existing themes. They represent pumps with entirely new capabilities (carotenoid antenna binding) and unparalleled functional plasticity (pump/sensor switching). As such, HeimdallRs and CryoRs are high-priority targets for structural biology to elucidate these novel mechanisms and fill the "Unknown" archaeal function category.

### 4.4. New Frontiers in Eukaryota: Filling the "Unknown" Gap

The eukaryotic domain is also yielding a wealth of novel opsin families, many with unknown or poorly understood functions.
**Apusomonad Rhodopsins (ApuRs):** A groundbreaking discovery from 2025, ApuRs are a completely new family of anion-selective channels found in apusomonads, a group of heterotrophic flagellate protists sister to the opisthokonts (which include animals and fungi).⁵⁰ The discovery is highly significant for two reasons. First, phylogenetic analysis shows that ApuRs evolved their channel function independently of all other known channelrhodopsin families, representing a remarkable case of convergent evolution for this complex function. Second, they possess the most blue-shifted absorption spectra of any known rhodopsin channel, with absorption maxima in the near-UV and violet range, suggesting adaptation to a unique light environment.⁵¹
**Xenopsins and Other Orphans:**

- **Xenopsins:** This is a phylogenetically distinct and poorly characterized opsin family found in protostomes, including mollusks and annelids.⁵⁹ While they are presumed to be involved in photoreception, their specific role is unclear. In a striking finding, a xenopsin in the chiton *Leptochiton asellus* was found to be co-expressed with a canonical rhabdomeric opsin (r-opsin) in the same larval photoreceptor cells, which themselves bear both microvilli and cilia.⁶⁰ This co-expression of two distinct opsin types hints at complex, unknown signaling or processing roles within a single photoreceptor.

- **Orphan Fungal Opsins:** While some fungal rhodopsins are now known to be proton pumps (e.g., from *Leptosphaeria maculans* and *Fusarium fujikuroi*)⁹, many others remain functionally uncharacterized. The opsin from the model organism *Neurospora crassa* (Neurospora Opsin-1, or NR), for example, has been known for over two decades, yet its physiological function remains enigmatic, though it has been implicated in regulating sexual reproduction.⁹ Genomic studies of fungi like *Neurospora* reveal hundreds of lineage-specific "orphan genes" with no known homologs or functions, representing a deep reservoir for the potential discovery of novel opsin families.⁶¹
For the Eukaryota "Unknown" category, ApuRs provide a spectacular example of a novel, structurally uncharacterized family with a known (but convergently evolved) function. Xenopsins, along with the vast number of orphan fungal and animal opsins identified in genomic surveys⁶⁷, represent a frontier of proteins with truly unknown physiological roles that are ripe for functional and structural investigation.

## V. Synthesis, Recommendations, and Future Outlook

### 5.1. A Revised Opsin Landscape: From Gaps to Hypotheses

This comprehensive literature review, prompted by the initial structural analysis, transforms the view of the microbial opsin landscape from a simple matrix with empty cells into a dynamic map shaped by clear evolutionary principles. The eight "missing" function-domain combinations are now resolved into either fillable data gaps (the archaeal chloride pump) or true biological absences governed by fundamental evolutionary divergences (prokaryotic channels and non-bacterial sodium pumps).
Several major themes emerge that explain this revised landscape. First is the fundamental pump-channel dichotomy, rooted in the opposing biophysical requirements for energy conservation versus rapid signaling, which drove the evolution of distinct mechanisms in prokaryotes (two-component sensory-pumps) and eukaryotes (single-component channels). Second is the pervasive role of Lateral Gene Transfer (LGT), which, coupled with gene loss, creates the "patchy" distribution of functions like the archaeal chloride pump, explaining its absence in limited datasets. Third is the power of convergent evolution, where the robust 7TM scaffold has been independently sculpted multiple times to achieve the same function (e.g., inward chloride pumping) using different molecular solutions.
These high-level principles are ultimately encoded at the molecular level. For ion pumps, function is largely determined by a small number of critical amino acid residues in the third transmembrane helix, forming a "motif" that dictates ion specificity. A synthesis of the literature provides a clear functional key for interpreting new sequences.
### Table 3: Key Amino Acid Motifs Determining Ion Specificity in Microbial Rhodopsin Pumps

| Function | Key Motif | Residue Positions (BR) | Domain(s) | Example Protein | References |
|----------|-----------|------------------------|-----------|-----------------|------------|
| Outward H⁺ Pump | DTD / DTE | 85, 89, 96 | Archaea, Bacteria, Eukarya | Bacteriorhodopsin (BR), Proteorhodopsin (PR) | 9 |
| Inward Cl⁻ Pump | TSA | 85, 89, 96 | Archaea | Halorhodopsin (HR) | 15 |
| Inward Cl⁻ Pump | NTQ | 85, 89, 96 | Bacteria | N. marinus Rhodopsin 3 (NM-R3) | 12 |
| Inward Cl⁻ Pump | TSD | 85, 89, 96 | Bacteria (Cyanobacteria) | Mastigocladopsis repens Rhodopsin (MastR) | 11 |
| Outward Na⁺ Pump | NDQ | 85, 89, 96 | Bacteria | Krokinobacter rhodopsin 2 (KR2) | 12 |
| Inward H⁺ Pump | FTD (example) | 85, 89, 96 | Archaea, Bacteria | Schizorhodopsin (SzR) | 23 |


### 5.2. Summary of Counterion Systems in Microbial Rhodopsins
Rhodopsin Class	Key Motif / Residues	Counterion System
Outward Proton Pumps		
Bacteriorhodopsin (BR)	DTD (Asp-Thr-Asp)	Complex of two Aspartate residues (Asp85 and Asp212).
Proteorhodopsin (PR)	DTE (Asp-Thr-Glu)	Complex of Aspartate and Glutamate residues (Asp85 and Glu96, BR numbering).
Inward Chloride Pumps		
Halorhodopsin (HR)	TSA (Thr-Ser-Ala)	The transported Chloride (Cl⁻) ion itself binds near the PRSB and serves as the counterion.
Bacterial Cl⁻ Pumps (NTQ)	NTQ (Asn-Thr-Gln)	The transported Chloride (Cl⁻) ion serves as the counterion.
Cyanobacterial Cl⁻ Pumps (TSD)	TSD (Thr-Ser-Asp)	The transported Chloride (Cl⁻) ion serves as the counterion.
Outward Sodium Pump
KR2	NDQ (Asn-Asp-Gln)	A single Aspartate residue (Asp116) from the NDQ motif acts as the counterion.
Channelrhodopsins		
Cation Channelrhodopsins (CCRs)	Conserved carboxylates	A complex of two carboxylate residues (Asp/Glu) and a surrounding network of Glutamates creates a negatively charged pore.
Anion Channelrhodopsins (ACRs)	Non-carboxylate at primary site	Lacks a primary carboxylate counterion. The protonated Schiff base itself is thought to directly mediate anion flux.
Sensory Rhodopsins		
SR-I & SR-II	Conserved Asp residues	A complex of two Aspartate residues (Asp75 and Asp201 in NpSRII), modulated by a nearby Arginine (Arg72).
Novel Families		
Heliorhodopsins (HeRs)	ESL (Glu-Ser-Leu)	A single Glutamate residue (e.g., Glu107 in HeR-48C12).
Kalium Channelrhodopsins (KCRs)	Asp105, Asp116	Asp105 and Asp116 are identified as key residues for K⁺ conductance, but are not explicitly defined as the counterion complex.
Apusomonad Rhodopsins (ApuRs)	DXQ / XTQ	The counterion system is not fully defined but involves divergent motifs. Some members have an Aspartate at the primary position, which is unusual for an anion channel.
CryoRhodopsins (CryoRs)	RRESEDK / RREAEDK	The specific counterion system is not detailed in the available literature.


### 5.2. A Roadmap for Structural and Functional Characterization

The clarified landscape points toward clear, actionable priorities for future experimental and computational work aimed at completing the structural and functional map of microbial opsins.
**High-Priority Targets for Structural Biology:**

1. **Heimdallarchaeial Rhodopsins (HeimdallRs):** Solving the high-resolution structure of a HeimdallR is of paramount importance. A structure, ideally captured in complex with its bound carotenoid antenna, would provide the first atomic-level view of an opsin from the critical Asgard archaeal lineage. It would reveal the molecular basis of the carotenoid-binding fenestration and the mechanism of energy transfer, offering profound insights into the evolution of light-harvesting and the metabolic capabilities of the archaeal ancestor of all eukaryotes.

2. **Apusomonad Rhodopsins (ApuRs):** Determining the structure of an ApuR would be a landmark achievement, revealing a third, independently evolved molecular solution for creating a light-gated ion channel. A comparative structural analysis of ApuR with chlorophyte CCRs and cryptophyte ACRs would be invaluable for uncovering the conserved and divergent principles of channel architecture, ion selectivity, and light-dependent gating.

3. **CryoRhodopsins (CryoRs):** Capturing the structures of CryoR in both its pump-competent and sensor-competent states would provide an unprecedented view of functional switching within a single protein scaffold. This would elucidate how subtle conformational changes, potentially modulated by the unique buried arginine residue, can toggle the protein between two distinct functional outputs.
**A Strategy for Navigating the "Unknown":**

The most effective path forward for characterizing the vast "Unknown" functional space is to systematically mine the ever-growing genomic and metagenomic databases for Alt-Rhodopsin (AltR) sequences. By definition, these proteins deviate from the canonical functional motifs. A screening strategy should prioritize sequences with novel substitutions at the critical Arg82 position, particularly those originating from under-sampled phylogenetic clades. These candidates should then be subjected to heterologous expression and high-throughput functional screening (e.g., pH or ion-sensitive fluorescence assays). The most promising candidates—those exhibiting novel ion specificities or other unexpected activities—should then be prioritized for detailed biophysical and structural characterization.

### 5.3. Concluding Perspective: The Expanding Opsin Universe

The integration of predictive structural biology with a deep literature review demonstrates that while our knowledge of the microbial opsin family has grown exponentially, we are far from a complete census. The landscape is not static; the rate of discovery of novel families with unexpected functions and evolutionary histories is accelerating, driven by the power of large-scale sequencing.⁵
The future of the field lies in a continued, tight integration of genomics, computational structure prediction, and high-throughput functional screening to navigate the enormous, uncharacterized sequence space. The most profound discoveries will likely not be mere variations on known themes, but entirely new functional paradigms, such as the light-dependent enzymatic regulation performed by Heliorhodopsins, which were hidden until recently within the genomic "dark matter" of the microbial world.³ The quest to complete the opsin structural map is not simply about filling in the cells of a pre-defined matrix, but about continually redrawing the boundaries of what we understand a rhodopsin can be and do.
## References

1. The Microbial Opsin Family of Optogenetic Tools - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC4166436/
2. Structure-based insights into evolution of rhodopsins - ResearchGate, Accessed June 24, 2025, https://www.researchgate.net/publication/352865633_Structure-based_insights_into_evolution_of_rhodopsins
3. The Evolutionary Kaleidoscope of Rhodopsins - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC9601134/
4. The Evolutionary Kaleidoscope of Rhodopsins | mSystems - ASM Journals, Accessed June 24, 2025, https://journals.asm.org/doi/10.1128/msystems.00405-22
5. Transporter Classification Database (TCDB): 2021 update | Nucleic Acids Research | Oxford Academic, Accessed June 24, 2025, https://academic.oup.com/nar/article/49/D1/D461/5973435
6. CryoRhodopsins: a comprehensive characterization of a group of..., Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2024.01.15.575777v2
7. Microbial Rhodopsins: Diversity, Mechanisms, and Optogenetic Applications - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC5747503/
8. Microbial rhodopsins: wide distribution, rich diversity and great potential - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC4736836/
9. Rhodopsins at a glance | Journal of Cell Science | The Company of Biologists, Accessed June 24, 2025, https://journals.biologists.com/jcs/article/134/22/jcs258989/273540/Rhodopsins-at-a-glance
10. Microbial Halorhodopsins: Light-Driven Chloride Pumps | Chemical Reviews, Accessed June 24, 2025, https://pubs.acs.org/doi/abs/10.1021/acs.chemrev.7b00715
11. The crystal structures of a chloride-pumping microbial rhodopsin and its proton-pumping mutant illuminate proton transfer determinants - PubMed Central, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC7606686/
12. Functional characterization of flavobacteria rhodopsins reveals a unique class of light-driven chloride pump in bacteria | PNAS, Accessed June 24, 2025, https://www.pnas.org/doi/10.1073/pnas.1403051111
13. Evolution of rhodopsin ion pumps in haloarchaea - PMC - PubMed Central, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC1885257/
14. Microbial rhodopsins: functional versatility and genetic mobility - PubMed, Accessed June 24, 2025, https://pubmed.ncbi.nlm.nih.gov/17008099/
15. Ion-pumping microbial rhodopsins - Frontiers, Accessed June 24, 2025, https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2015.00052/full
16. A, positions of three amino acid residues consisting of the motif.... - ResearchGate, Accessed June 24, 2025, https://www.researchgate.net/figure/A-positions-of-three-amino-acid-residues-consisting-of-the-motif-X-ray-crystal_fig2_284161438
17. Functional characterization of flavobacteria rhodopsins reveals a unique class of light-driven chloride pump in bacteria | Request PDF - ResearchGate, Accessed June 24, 2025, https://www.researchgate.net/publication/261409809_Functional_characterization_of_flavobacteria_rhodopsins_reveals_a_unique_class_of_light-driven_chloride_pump_in_bacteria
18. Functional characterization of flavobacteria rhodopsins reveals a unique class of light-driven chloride pump in bacteria - DSpace@MIT, Accessed June 24, 2025, https://dspace.mit.edu/handle/1721.1/91985?show=full
19. Ion-pumping microbial rhodopsins - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC4585134/
20. pmc.ncbi.nlm.nih.gov, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC6457933/#:~:text=The%20first%20microbial%20rhodopsin%20(bacteriorhodopsin,for%20the%20cell%20(2).
21. Energetics and dynamics of a light-driven sodium-pumping rhodopsin - PNAS, Accessed June 24, 2025, https://www.pnas.org/doi/10.1073/pnas.1703625114
22. Discovery of a new light-powered sodium pump | The University of Tokyo, Accessed June 24, 2025, https://www.u-tokyo.ac.jp/focus/en/articles/a_00197.html
23. Ion-transporting mechanism in microbial rhodopsins: Mini-review relating to the session 5 at the 19th International Conference on Retinal Proteins, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC10865854/
24. Structure and mechanisms of sodium-pumping KR2 rhodopsin - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC6457933/
25. Na + Binding and Transport: Insights from Light-Driven Na + -Pumping Rhodopsin - MDPI, Accessed June 24, 2025, https://www.mdpi.com/1420-3049/28/20/7135
26. Ion-transporting mechanism in microbial rhodopsins: Mini-review..., Accessed June 24, 2025, https://www.jstage.jst.go.jp/article/biophysico/20/Supplemental/20_e201005/_html/-char/en
27. 1 Microbial Rhodopsins: Phylogenetic and Functional Diversity - CiteSeerX, Accessed June 24, 2025, https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=129bb9335f53a70e8afe07d5b91de6226893f14a
28. Light-harvesting by antenna-containing rhodopsins in pelagic Asgard archaea - bioRxiv, Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2024.09.18.613612v1.full.pdf
29. Structural insights into light harvesting by antenna-containing..., Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC12137139/
30. Proton-pumping rhodopsins in marine diatoms - bioRxiv, Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2022.01.18.476826v1.full.pdf
31. Leptosphaeria rhodopsin: Bacteriorhodopsin-like proton pump from a eukaryote | PNAS, Accessed June 24, 2025, https://www.pnas.org/doi/10.1073/pnas.0409659102
32. Ion channels versus ion pumps: the principal difference, in principle..., Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC2742554/
33. The Origin and Early Evolution of Membrane Channels - ResearchGate, Accessed June 24, 2025, https://www.researchgate.net/publication/8022340_The_Origin_and_Early_Evolution_of_Membrane_Channels
34. Optogenetics Guide - Addgene, Accessed June 24, 2025, https://www.addgene.org/guides/optogenetics/
35. Channelrhodopsins: light-activated ion channels - — Institut für Biologie - Humboldt-Universität zu Berlin, Accessed June 24, 2025, https://www.biologie.hu-berlin.de/de/gruppenseiten/expbp/chanelrhodopsin
36. Emerging Diversity of Channelrhodopsins and Their Structure-Function Relationships - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC8818676/
37. Cation and Anion Channelrhodopsins: Sequence Motifs and Taxonomic Distribution - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC8406140/
38. Channelrhodopsin - Wikipedia, Accessed June 24, 2025, https://en.wikipedia.org/wiki/Channelrhodopsin
39. Microbial and Animal Rhodopsins: Structures, Functions, and Molecular Mechanisms | Chemical Reviews - ACS Publications, Accessed June 24, 2025, https://pubs.acs.org/doi/10.1021/cr4003769
40. Channelrhodopsin-2, a directly light-gated cation-selective membrane channel | PNAS, Accessed June 24, 2025, https://www.pnas.org/doi/abs/10.1073/pnas.1936192100?doi=10.1073/pnas.1936192100
41. Channelrhodopsin unchained: Structure and mechanism of a light-gated cation channel - Refubium - Freie Universität Berlin, Accessed June 24, 2025, https://refubium.fu-berlin.de/bitstream/fub188/14801/1/BBABIO-13-127R1_MS.pdf
42. Channelrhodopsin-2, a directly light-gated cation-selective membrane channel - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC283525/
43. Sodium-Selective Channelrhodopsins - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC11592924/
44. Robust Optogenetic Inhibition with Red-light-sensitive Anion-conducting Channelrhodopsins - eLife, Accessed June 24, 2025, https://elifesciences.org/reviewed-preprints/90100v1
45. Structural Changes in an Anion Channelrhodopsin: Formation of the K and L Intermediates at 80 K, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC5747504/
46. RubyACRs, nonalgal anion channelrhodopsins with highly red-shifted absorption | PNAS, Accessed June 24, 2025, https://www.pnas.org/doi/10.1073/pnas.2005981117
47. Gating mechanisms of a natural anion channelrhodopsin - PNAS, Accessed June 24, 2025, https://www.pnas.org/doi/10.1073/pnas.1513602112
48. Anion-conducting channelrhodopsin - Wikipedia, Accessed June 24, 2025, https://en.wikipedia.org/wiki/Anion-conducting_channelrhodopsin
49. A phylogenetic tree of three families of channelrhodopsins constructed... - ResearchGate, Accessed June 24, 2025, https://www.researchgate.net/figure/A-phylogenetic-tree-of-three-families-of-channelrhodopsins-constructed-by-the_fig1_320651153
50. Apusomonad rhodopsins, a new family of ultraviolet to blue light absorbing rhodopsin channels - bioRxiv, Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2025.04.02.646541v1.full.pdf
51. Apusomonad rhodopsins, a new family of ultraviolet to blue light..., Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2025.04.02.646541v1
52. Light at the end of the tunnel of Corti - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC10238449/
53. The Ion-translocating Microbial Rhodopsin (MR) Family - TCDB » SEARCH, Accessed June 24, 2025, https://tcdb.org/search/result.php?tc=3.E.1
54. Proton- transporting heliorhodopsins from marine giant viruses - eLife, Accessed June 24, 2025, https://elifesciences.org/articles/78416.pdf
55. No heliorhodopsins in Gram-negative bacteria | Research Communities by Springer Nature, Accessed June 24, 2025, https://communities.springernature.com/posts/no-heliorhodopsins-in-gram-negative-bacteria
56. Heliorhodopsin Helps Photolyase to Enhance the DNA Repair Capacity | Microbiology Spectrum - ASM Journals, Accessed June 24, 2025, https://journals.asm.org/doi/10.1128/spectrum.02215-22
57. Heliorhodopsin binds and regulates glutamine synthetase activity..., Accessed June 24, 2025, https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001817
58. Light-harvesting by antenna-containing rhodopsins in pelagic Asgard archaea | bioRxiv, Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2024.09.18.613612v1.full-text
59. Molluscan Genomes Reveal Extensive Differences in Photopigment Evolution Across the Phylum - Oxford Academic, Accessed June 24, 2025, https://academic.oup.com/mbe/article/40/12/msad263/7457512
60. Co-expression of xenopsin and rhabdomeric opsin in photoreceptors bearing microvilli and cilia | eLife, Accessed June 24, 2025, https://elifesciences.org/articles/23435
61. Orphan genes are clustered with allorecognition loci and may be involved in incompatibility and speciation in Neurospora | bioRxiv, Accessed June 24, 2025, https://www.biorxiv.org/content/10.1101/2022.06.10.495464v1.full
62. Heliorhodopsin Evolution Is Driven by Photosensory Promiscuity in Monoderms | mSphere, Accessed June 24, 2025, https://journals.asm.org/doi/10.1128/mSphere.00661-21
63. Heliorhodopsin - Wikipedia, Accessed June 24, 2025, https://en.wikipedia.org/wiki/Heliorhodopsin
64. High-resolution structural insights into the heliorhodopsin family - PMC, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC7049168/
65. The Evolutionary Kaleidoscope of Rhodopsins | mSystems - ASM Journals, Accessed June 24, 2025, https://journals.asm.org/doi/abs/10.1128/msystems.00405-22
66. Orphan genes are clustered with allorecognition loci and may be involved in incompatibility and speciation in Neurospora - ResearchGate, Accessed June 24, 2025, https://www.researchgate.net/publication/361261175_Orphan_genes_are_clustered_with_allorecognition_loci_and_may_be_involved_in_incompatibility_and_speciation_in_Neurospora
67. Shedding new light on opsin evolution - PMC - PubMed Central, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC3223661/
68. Opsin expression varies across larval development and taxa in pteriomorphian bivalves, Accessed June 24, 2025, https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2024.1357873/full
69. Orphan GPCRs and neuromodulation - PMC - PubMed Central, Accessed June 24, 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC3474844/
70. Taxonomic distribution of opsin families inferred from UniProt Reference Proteomes and a suite of opsin-specific hidden Markov models - Frontiers, Accessed June 24, 2025, https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2023.1190549/full
