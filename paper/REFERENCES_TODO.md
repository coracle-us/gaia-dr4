# References for known_companion_predictions.tex

DONE 2026-11-24: known_companion_predictions.bib is now the ADS BibTeX export (fetched with the ADS API
for the 10 bibcodes, then rekeyed by scripts/ads_rekey_bib.py). To refresh it, repeat the steps below:

1. Paste this into the ADS search box (ui.adsabs.harvard.edu):

   identifier:("arXiv:1312.2954" OR "arXiv:1601.07595" OR "arXiv:2309.08665" OR "arXiv:2511.04673" OR "arXiv:2609.09077" OR "arXiv:2610.25152" OR "arXiv:2609.01377" OR "arXiv:2206.05439" OR "arXiv:2610.07149" OR "arXiv:2105.11662")

   Expect 10 records. If fewer come back, search the missing arXiv ids one at a time.
2. Select all, Export -> BibTeX, save the file.
3. python scripts/ads_rekey_bib.py <saved file>
4. Upload paper/known_companion_predictions.tex, .bib and figures/known_companion_predictions_fig1.pdf
   to Overleaf (AASTeX 7 template, aastex701.cls + aasjournalv7.bst) and compile.

| Key | arXiv | Published as (check against ADS) |
|---|---|---|
| knutson2014 | 1312.2954 | ApJ 785, 126 |
| bryan2016 | 1601.07595 | ApJ 821, 89 |
| espinoza2023 | 2309.08665 | AJ 166, 231 |
| lammers2025 | 2511.04673 | AJ 171, 18 (2026; key kept) |
| ranshaw2026 | 2609.09077 | AJ 172, 254 |
| parc2026 | 2610.25152 | preprint |
| xiao2026 | 2609.01377 | ApJL 1009, L6 |
| holl2023 | 2206.05439 | A&A 674, A10 |
| greenfield2026 | 2610.07149 | preprint |
| brandt2021 | 2105.11662 | ApJS 254, 42 |

Nothing in this note needs the ESA Gaia archive: the host parameters come from data/targets.csv
(all 36 hosts complete) and the scanning law from the local gaiascanlaw data.
