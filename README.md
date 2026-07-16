<p align="center">
  <img src="img/esKaKs.svg" height="200" alt="eskaks logo" />
</p>

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-%23af64d1?style=flat-square)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-%23dea584?style=flat-square)](https://www.rust-lang.org/)
[![Version](https://img.shields.io/badge/version-1.4.0-%23149389?style=flat-square)](https://github.com/PathoGenOmics-Lab/eskaks/releases)
[![CI](https://github.com/PathoGenOmics-Lab/eskaks/actions/workflows/ci.yml/badge.svg)](https://github.com/PathoGenOmics-Lab/eskaks/actions/workflows/ci.yml)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square)](https://github.com/PathoGenOmics-Lab)

**Pairwise dN/dS (Ka/Ks) and per-gene pN/pS from codon-aligned sequences or VCF files.**
**Pure Rust · no C dependencies · self-contained interactive HTML reports.**

[Tutorial](docs/tutorial.md) · [Get started](#get-started) · [Docs](docs/) · [Citation](#citation)

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>

---

## What is eskaks?

Fast **pairwise dN/dS (Ka/Ks)** and **per-gene pN/pS** for measuring natural selection
on protein-coding genes. Two modes:

- **`eskaks fasta`**: pairwise dN/dS from codon-aligned FASTA (Nei-Gojobori or Li/LPB93 model).
- **`eskaks vcf`**: per-gene pN/pS from a VCF + reference + GFF3, with a genome-wide
  neutrality scan and a self-contained interactive HTML report.

Pure Rust, no C dependencies, up to **2,641× faster** than KaKs_Calculator (R² = 1.0
on the Li model). New here? Start with the **[hands-on tutorial](docs/tutorial.md)**.

## Get started

Requires [Rust](https://www.rust-lang.org/tools/install) ≥ 1.70.

```bash
git clone https://github.com/PathoGenOmics-Lab/eskaks.git
cd eskaks && make release && cp target/release/eskaks ~/.local/bin/
```

Try it on the bundled [`examples/`](examples/), no data of your own needed:

```bash
# pairwise dN/dS from aligned sequences
eskaks fasta examples/genes.fasta -o first_run

# per-gene pN/pS + an interactive HTML report
eskaks vcf --ref examples/toy_genome/reference.fasta \
  --gff examples/toy_genome/genes.gff3 --vcf examples/toy_genome/variants.vcf \
  --genetic-code 11 --report -o toy_scan   # → open toy_scan_report.html
```

`--report` builds a single self-contained HTML dashboard (interactive Manhattan /
volcano / QQ, McDonald-Kreitman, colour-blind mode, scales to whole genomes; no
internet needed). The [tutorial](docs/tutorial.md) walks through reading the output.

## Documentation

Full docs live in **[`docs/`](docs/)**: [tutorial](docs/tutorial.md) ·
[CLI reference](docs/cli-reference.md) · [VCF analysis](docs/vcf-analysis.md) ·
[interpreting results](docs/interpreting-results.md) · [glossary](docs/glossary.md) ·
[performance](docs/performance.md) · [FAQ](docs/faq.md).

## Citation

If you use eskaks in your research, please cite:

> Ruiz-Rodriguez P, Coscollá M. **eskaks: fast pairwise dN/dS and per-gene pN/pS from sequences or VCFs.** https://github.com/PathoGenOmics-Lab/eskaks

```bibtex
@software{ruiz-rodriguez_eskaks_2026,
  title   = {eskaks: fast pairwise dN/dS and per-gene pN/pS from sequences or VCFs},
  author  = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  year    = {2026},
  url     = {https://github.com/PathoGenOmics-Lab/eskaks},
  version = {1.4.0},
  license = {GPL-3.0}
}
```

## License

[GNU General Public License v3.0](LICENSE)

---

<h2 id="contributors" align="center">

✨ Contributors
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
eskaks is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Design">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
