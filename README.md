# annotations
Repository for VEP annotations and "extra" annotations.

Each ```v0.x``` is associated with a published Dockstore workflow. 

**```vepAnnotateHail_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/annotations/vep-annotate-hail-v01
- Relatively stable. Currently uses VEP v112 but **TODO**: edit so using different VEP version inputs only involves changing a few WDL inputs.
- Plugins:
  - AlphaMissense
  - EVE 
  
**```vepAnnotateHailMT_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/annotations/vep-annotate-hail-mt-v01
- Relatively stable. Basically replicates ```vepAnnotateHail_v0.1.wdl``` but expects MT input and saves as MT (**TODO**: maybe output as VCF?).
- Rarely used, only if we have input in MT form and not VCF(s).

**```vepAnnotateHailExtra_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/annotations/vep-annotate-hail-extra-v01
- Potentially misleading name because VEP is not run, but this workflow is intended to be run after **```vep-annotate-hail-v(x.x)```** for "extra" annotations:
  - [Optional] Noncoding (from bed file)
  - LOEUF from gnomAD v2 and gnomAD v4
  - OMIM
  - REVEL
  - MPC
  - ClinVar
  - [Optional] Gene list(s)
  - SpliceAI
- ```v0.1```: annotate with only one (optional) gene list.
- ```v0.2```: annotate with any number of (optional) gene lists.
