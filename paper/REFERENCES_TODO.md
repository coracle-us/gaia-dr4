# References for known_companion_predictions.tex

The .bib in this folder is a placeholder. Replace it with ADS exports before compiling:

1. Paste this into the ADS search box (ui.adsabs.harvard.edu):

   identifier:("arXiv:1312.2954" OR "arXiv:1601.07595" OR "arXiv:2309.08665" OR "arXiv:2511.04673" OR "arXiv:2609.09077" OR "arXiv:2610.25152" OR "arXiv:2609.01377" OR "arXiv:2206.05439" OR "arXiv:2610.07149")

   Expect 9 records. If fewer come back, search the missing arXiv ids one at a time.
2. Select all, Export -> BibTeX, save the file.
3. python scripts/ads_rekey_bib.py <saved file>
4. Upload paper/known_companion_predictions.tex, .bib and figures/known_companion_predictions.pdf
   to Overleaf (AASTeX 7 template, aastex701.cls + aasjournalv7.bst) and compile.

| Key | arXiv | Published as (check against ADS) |
|---|---|---|
| knutson2014 | 1312.2954 | ApJ 785, 126 |
| bryan2016 | 1601.07595 | ApJ 821, 89 |
| espinoza2023 | 2309.08665 | AJ 166, 231 |
| lammers2025 | 2511.04673 | AJ, accepted |
| ranshaw2026 | 2609.09077 | AJ, accepted |
| parc2026 | 2610.25152 | preprint |
| xiao2026 | 2609.01377 | ApJL, accepted |
| holl2023 | 2206.05439 | A&A 674, A10 |
| greenfield2026 | 2610.07149 | preprint |

Nothing in this note needs the ESA Gaia archive: the host parameters come from data/targets.csv
(all 36 hosts complete) and the scanning law from the local gaiascanlaw data.
