
from enum import Enum
from pathlib import Path
from latch import large_task, workflow
from natsort import as_ascii

from nupack import *  # Import NUPACK

import matplotlib.pyplot as plt
import pandas as pd
import re
from latch.types import LatchFile, LatchMetadata, LatchAuthor, LatchParameter, LatchAppearanceType, LatchRule


class Material(Enum):
    dna = "DNA"
    rna = "RNA"
    rna95 = "RNA95"


class Ensemble(Enum):
    stacking = "stacking"
    nostacking = "nostacking"


@large_task
def tubeDesign(
    domains_file: LatchFile,
    target_strands_file: LatchFile,
    target_complexes_file: LatchFile,
    max_size: int = 3,
    tube_name: str = "t1",
    material: str = Material.rna,
    ensemble: str = Ensemble.stacking,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    out: str = "nupack-design",
) -> LatchFile:

    model = Model(material=material, ensemble=ensemble, celsius=temperature, sodium=sodium, magnesium=magnesium)

    degenerate_codes = r'[ATGCUMRWSYKVHDBN\d]+'
    splitter = r'(\S+)'
    concentration_pattern = r'>\s?(.*)\s+([0-9]{1,}\.?([eE][+-]|\.)[0-9]{1,})'
    dot_bracket_pattern = r'(^[\(\)\.+\d?\s?D?U?]+)'

    with open(domains_file, "r") as f:
        domains = f.read().splitlines()
    f.close()
    
    with open(target_strands_file, "r") as f:
        targetStrands = f.read().splitlines()
    f.close()

    with open(target_complexes_file, "r") as f:
        targetComplexes = f.read().splitlines()
    f.close()

    domainNames = [l[1:].strip() for l in domains if l.startswith('>')]
    domainSequences = [l.strip() for l in domains if re.match(degenerate_codes, l) and l.startswith('>') == False]
    domainDict = dict(zip(domainNames,domainSequences))

    targetStrandNames = [l[1:].strip() for l in targetStrands if l.startswith('>')]
    targetStrandDomains = [re.findall(splitter,l) for l in targetStrands if l.startswith('>') == False]
    targetStrandDomains = [[Domain(domainDict[el], name=el) for el in li] for li in targetStrandDomains]
    targetStrandDict = dict(zip(targetStrandNames,targetStrandDomains))

    targetComplexDetails = [
    {
        'complex': re.match(concentration_pattern,l).group(1),
        'conc': float(re.match(concentration_pattern,l).group(2)),
        'strands': re.findall(splitter,targetComplexes[targetComplexes.index(l)+1]),
        'structure' : re.match(dot_bracket_pattern, targetComplexes[targetComplexes.index(l)+2]).group(1),
    }
    for l in targetComplexes
    if l.startswith('>')]

    for l in range(len(targetComplexDetails)):
        targetComplexDetails[l]['strands'] = [TargetStrand(targetStrandDict[el],name=el) for el in targetComplexDetails[l]['strands']]
        targetComplexDetails[l]['complex'] = TargetComplex(targetComplexDetails[l]['strands'], targetComplexDetails[l]['structure'], name=targetComplexDetails[l]['complex'])

    targetComplexDict = {targetComplexDetails[l]['complex']:targetComplexDetails[l]['conc'] for l in range(len(targetComplexDetails))}
    targetTube = TargetTube(on_targets=targetComplexDict, name=tube_name, off_targets=SetSpec(max_size=max_size))

    design = tube_design(tubes=[targetTube], hard_constraints=[], soft_constraints=[], defect_weights=None, options=None, model=model)

    result = design.run(trials=1)[0]

    outFile = f"/{out}.txt"

    content = f"""
    NUPACK DESIGN RESULTS
    {result}
    """
    with open(outFile, "w") as f:
        f.write(content)

    return LatchFile(outFile, f"latch:///{outFile}")


metadata = LatchMetadata(
    display_name="NUPACK - Tube Design",
    documentation="https://docs.nupack.org/design",
    author=LatchAuthor(
        name="Shivaramakrishna Srinivasan",
        email="shivaramakrishna.srinivasan@gmail.com",
        github="https://github.com/shivaramakrishna99",
    ),
    repository="https://github.com/beliveau-lab/NUPACK",
    license="BSD-3-Clause",
)

metadata.parameters["domains_file"] = LatchParameter(
    display_name="List of domains as FASTA",
    description="File containing list of domain names and sequences in FASTA format",
    section_title="Input",
)
metadata.parameters["target_strands_file"] = LatchParameter(
    display_name="List of target strands as FASTA",
    description="File containing list of target strand names with their constituent in FASTA format",
)
metadata.parameters["target_complexes_file"] = LatchParameter(
    display_name="List of target complexes as FASTA",
    description="File containing list of target complex names with interacting strands and target secondary structure in FASTA format",
)
metadata.parameters["max_size"] = LatchParameter(
    display_name="Maximum Complex Size",
    description="Specify maximum number of interactions between strands",
    section_title="Tube Specifications",
)
metadata.parameters["tube_name"] = LatchParameter(
    display_name="Name of Tube",
    description="Provide a name for the tube"
)
metadata.parameters["material"] = LatchParameter(
    display_name="Nucleic Acid Type",
    description="Choose between DNA and RNA free energy parameter sets. Default is 'rna', based on Matthews et al., 1999",
    section_title="Model Specification",
    hidden=True,
)
metadata.parameters["ensemble"] = LatchParameter(
    display_name="Ensemble Type",
    description="Choose between stacking and non stacking ensemble states. Default is 'stacking'",
    hidden=True
)
metadata.parameters["temperature"] = LatchParameter(
    display_name="Temperature (in °C)",
    description="Temperature of system. Default: 37.0",
    hidden=True,
)
metadata.parameters["sodium"] = LatchParameter(
    display_name="Na+ (in M)",
    description="The total concentration of (monovalent) sodium, potassium, and ammonium ions, specified as molarity. Default: 1.0, Range: [0.05,1.1]",
    hidden=True,
    section_title="Additional Model Specification"
)
metadata.parameters["magnesium"] = LatchParameter(
    display_name="Mg++ (in nM)",
    description="The total concentration of (divalent) magnesium ions, specified as molarity. Default: 0.0, Range: [0.0,0.2]",
    hidden=True
)
metadata.parameters["out"] = LatchParameter(
    display_name="Output File Name",
    section_title="Output"
)


@workflow(metadata)
def tubeDesignNUPACK(
    domains_file: LatchFile,
    target_strands_file: LatchFile,
    target_complexes_file: LatchFile,
    max_size: int = 3,
    tube_name: str = "t1",
    material: Material = Material.rna,
    ensemble: Ensemble = Ensemble.stacking,
    temperature: float = 37.0,
    sodium: float = 1.0,
    magnesium: float = 0.0,
    out: str = "tube-design",
) -> LatchFile:
    """Design a tube containing multiple nucleic acid strands

    # NUPACK - Tube Design
    ---

    ## **How to use**
    ---

    1. Provide a list of domains required for your design in FASTA format. Make sure the domain name is a single word (no spaces) and the sequence is provided as per the degenerate nucleotide code scheme
    ```
    >domain_a
    A4
    >domain_b
    R5N5
    ```

    2. Provide a list of target strands and the constituent domains that make them up, in FASTA format. Make sure the strand name is a single word (no spaces) and list of domain names for the strand are separated by a space
    ```
    >strand_a
    a a
    >domain_b
    b a b
    ```

    3. Provide a list of target complexes, and concentrations, their constituent strands, and the targeted secondary structure, in FASTA format. Make sure the strand name is a single word (no spaces), and separate by a space and specify the concentration value for the complex in floating point notation. In the next line, specify a list of strand names for the complex separated by a space. In the line after that, include the intented secondary structure in DU+ notation 
    ```
    >C1 1e-8
    A B C
    ........((((((((((+))))))))))((((((((((+))))))))))..............
    >C2 1e-8
    B C
    .10(10+)10.14
    ```

    4. Set a maximum complex size for off-targets and give your tube design job a name.

    5. Specify any other changes in the construction of the Model() object using the hidden parameters such as ensemble type and ion concentrations. 

    6. Run the workflow!

    ## **About**
    ---

    [NUPACK](https://docs.nupack.org/#about) is a growing software suite for the analysis and design of nucleic acid structures, devices, and systems serving the needs of researchers in the fields of nucleic acid nanotechnology, molecular programming, synthetic biology, and across the life sciences more broadly.

     ## **Citations**
    ---

    ### NUPACK Design Algorithms

    - **Multi-tube design**
        - B. R. Wolfe, N. J. Porubsky, J. N. Zadeh, R. M. Dirks, and N. A. Pierce. Constrained multistate sequence design for nucleic acid reaction pathway engineering. [*J Am Chem Soc*](http://pubs.acs.org/doi/abs/10.1021/jacs.6b12693), 139:3134-3144, 2017. ([pdf](http://www.nupack.org/downloads/serve_public_file/wolfe17.pdf?type=pdf), [supp info](http://www.nupack.org/downloads/serve_public_file/wolfe17_supp.pdf?type=pdf))

    - **Test tube design**
        - B. R. Wolfe and N. A. Pierce. Sequence design for a test tube of interacting nucleic acid strands. [*ACS Synth Biol*](http://pubs.acs.org/doi/abs/10.1021/sb5002196), 4:1086-1100, 2015. ([pdf](http://www.nupack.org/downloads/serve_public_file/sb5002196.pdf?type=pdf), [supp info](http://www.nupack.org/downloads/serve_public_file/sb5002196_si_001.pdf?type=pdf), [supp tests](http://www.nupack.org/downloads/serve_public_file/sb5002196_si_002.zip?type=zip))

    - **Complex design**
        - J. N. Zadeh, B. R. Wolfe, and N. A. Pierce. Nucleic acid sequence design via efficient ensemble defect optimization. [*J Comput Chem*](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21633/abstract), 32:439--452, 2011. ([pdf](http://www.nupack.org/downloads/serve_public_file/jcc11b.pdf?type=pdf), [supp info](http://www.nupack.org/downloads/serve_public_file/JCC_21633_sm_suppinfo.pdf?type=pdf), [supp tests](http://www.nupack.org/downloads/serve_public_file/JCC_21633_sm_supptests.zip?type=zip))

    - **Design paradigms**
        - R. M. Dirks, M. Lin, E. Winfree, and N. A. Pierce. Paradigms for computational nucleic acid design. [*Nucl Acids Res*](https://academic.oup.com/nar/article/32/4/1392/1038453), 32:1392-1403, 2004. ([pdf](http://www.nupack.org/downloads/serve_public_file/nar04.pdf?type=pdf), [supp info](http://www.nupack.org/downloads/serve_public_file/nar04_supp.pdf?type=pdf), [supp seqs](http://www.nupack.org/downloads/serve_public_file/nar04_sequences.tar.gz?type=zip))

    **Acknowledgements** - (https://docs.nupack.org/#acknowledgments)

    *Authored by Shivaramakrishna Srinivasan. Feel free to reach out to me at shivaramakrishna.srinivasan@gmail.com*
    ---
    """

    return tubeDesign(
        domains_file=domains_file,
        target_strands_file=target_strands_file,
        target_complexes_file=target_complexes_file,
        max_size=max_size,
        tube_name=tube_name,
        material=material,
        ensemble=ensemble,
        temperature=temperature,
        sodium=sodium,
        magnesium=magnesium,
        out=out
    )
