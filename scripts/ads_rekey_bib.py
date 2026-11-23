"""Turn an ADS BibTeX export into paper/known_companion_predictions.bib with the note's citation keys.

1. On ADS (ui.adsabs.harvard.edu) run the query in paper/REFERENCES_TODO.md, select all,
   Export -> BibTeX, save the file (e.g. ~/Downloads/export-bibtex.bib).
2. python scripts/ads_rekey_bib.py ~/Downloads/export-bibtex.bib
The script matches each ADS record to a citation key by its arXiv eprint, writes the bib, and
fails if any key cited in the .tex is missing. ADS entries are written unedited apart from the key.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "paper" / "known_companion_predictions.tex"
OUT = ROOT / "paper" / "known_companion_predictions.bib"
KEYS = {  # arXiv id -> key used in the .tex
    "1312.2954": "knutson2014", "1601.07595": "bryan2016", "2309.08665": "espinoza2023",
    "2511.04673": "lammers2025", "2609.09077": "ranshaw2026", "2610.25152": "parc2026",
    "2609.01377": "xiao2026", "2206.05439": "holl2023", "2610.07149": "greenfield2026",
}


def main(export):
    text = Path(export).read_text()
    entries = re.split(r"\n(?=@)", text.strip())
    out, found = [], set()
    for e in entries:
        m = re.search(r"eprint\s*=\s*\{(?:arXiv:)?([0-9]{4}\.[0-9]{4,5})", e)
        if not m or m.group(1) not in KEYS:
            print("skipped (no matching arXiv id):", e.splitlines()[0])
            continue
        key = KEYS[m.group(1)]
        out.append(re.sub(r"^(@\w+\{)[^,]+,", rf"\g<1>{key},", e, count=1))
        found.add(key)
    cited = set()
    for grp in re.findall(r"\\cite[a-z]*\{([^}]*)\}", TEX.read_text()):
        cited.update(k.strip() for k in grp.split(","))
    missing = sorted(cited - found)
    if missing:
        sys.exit(f"Not in the ADS export: {missing}. Search those arXiv ids on ADS and export again.")
    OUT.write_text("% ADS BibTeX export, rekeyed by scripts/ads_rekey_bib.py\n\n" + "\n\n".join(out) + "\n")
    print(f"wrote {OUT} with {len(out)} entries; all {len(cited)} cited keys present")


if __name__ == "__main__":
    main(sys.argv[1])
