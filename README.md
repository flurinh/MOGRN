<h1 align="center">MOGRN</h1>

<p align="center">
  <b>Microbial Opsin Generic Residue Numbering.</b><br>
  One standardized coordinate system across hundreds of diverse type-I opsins — so any two can be compared, residue for residue.
</p>

<p align="center"><img src="docs/grn-positions.jpg" alt="Key GRN microswitch positions on bacteriorhodopsin" width="380"></p>

<p align="center">
  <a href="https://flurinh.github.io/aboutme">◆ Portfolio</a> &nbsp;·&nbsp;
  <b>The build:</b>
  <a href="https://github.com/flurinh/LM-DTA">LM-DTA</a> →
  <a href="https://github.com/flurinh/mt">Master thesis</a> →
  <a href="https://github.com/flurinh/protos">ProtOS</a> →
  <b>MOGRN</b> →
  <a href="https://github.com/flurinh/lambda">Lambda</a> →
  <a href="https://github.com/flurinh/Protos_MCP">ProtOS-MCP</a>
</p>

---

## What it is

Microbial opsins are light-driven pumps, channels, and sensors built on a shared seven-helix
fold. MOGRN assigns a **Generic Residue Numbering** (GRN) that maps each residue to a
structurally equivalent position across the whole family, then uses it for **conservation
analysis, motif detection, co-evolution, and pump-vs-channel function discrimination** at
scale. It is the basis of my lead-author residue-numbering paper.

Anchor positions (the `.50` set): **3.50** — the functional switch (T in pumps → C in
channels); **6.50** — the conserved retinal-pocket tryptophan; **7.50** — the Schiff-base lysine.

<p align="center"><img src="docs/function-domain.png" alt="Molecular function by domain of life across the type-I opsin set" width="640"></p>
<p align="center"><i>Molecular function × domain of life across the type-I opsin set.</i></p>

## Built on ProtOS

MOGRN runs on the **[ProtOS](https://github.com/flurinh/protos)** framework (GRN processor,
structure/sequence handling):

```bash
git clone https://github.com/flurinh/protos.git && pip install -e protos
git clone https://github.com/flurinh/MOGRN.git && cd MOGRN
pip install -r requirements.txt
```

## Pipeline

```bash
python prepare_data.py    # 1. validate structures, set up caching
python prepare_yaml.py    # 2. sequence configuration
# 3–6: GRN assignment → conservation → motif detection → function analysis
```

See [`GUIDE.md`](GUIDE.md) for the full walkthrough and the required data folders.

## How it fits

GRN is the coordinate system that lets **[Lambda](https://github.com/flurinh/lambda)** line up
binding pockets across opsins and predict their colour.

---

<p align="center">
◀ <b>Previously:</b> <a href="https://github.com/flurinh/protos">ProtOS — the framework underneath</a>
&nbsp;·&nbsp;
<b>Next:</b> <a href="https://github.com/flurinh/lambda">Lambda — predicting opsin colour</a> ▶
</p>
