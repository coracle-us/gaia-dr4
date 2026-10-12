# AAS journals: AI policy and RNAAS rules (checked 2026-10-12)

## AI guidelines for authors

Source: AAS Publications Committee, "Guidelines on the use of AI / LLMs in preparing manuscripts" (authors), first posted **9 Sep 2026**; announcement https://journals.aas.org/?p=4629 ; full text https://journals.aas.org/author-llm-guidelines ; linked from https://journals.aas.org/policies/ .
- They apply to new submissions "later this fall"; no exact date is given.
- They are a living document. **Re-read them before submitting.**

| Topic | What the guidelines say |
|---|---|
| **Authorship** | AI cannot be an author (it cannot take responsibility for the content). |
| **Disclosure** | Required on every submission: state whether AI/LLMs were used, naming the tool and its purpose. It goes in an "AI disclosure statement included at the end of the article". Optionally, tell the editor which sections used it. |
| **Allowed, with disclosure and supervision** | Language editing, translation (checked by an author), accessibility, data/code formatting. Drafting text is allowed only if the authors read, check and edit all of it ("not generated unsupervised"). Explanatory illustrations: disclose the tool and check for copying. Data visualisations: check faithfulness to the data; the caption names the tool. |
| **Analysis code / data analysis** | The guidelines make no recommendation on AI use in data analysis or in writing analysis code. |
| **Never allowed** | AI content presented as data or as results from non-AI sources. Undisclosed AI use on figures, data, text or references. Plagiarism via AI. |
| **References** | AI must not be the sole source of the reference list. Check that every reference exists and is relevant (download it; use ADS). Fictitious references can lead to retraction. |
| **Responsibility** | The authors are fully accountable for all figures, data, text and references. |

**Consequences for this project**
1. Claude is not an author. Jon McBee is the sole author, so the Note uses "I" (RNAAS style for single-author Notes).
2. The Note ends with an AI disclosure statement. A draft is in `paper/dr4_companions.tex`. It must stay true:
   - Claude wrote most of the analysis code and drafts of text.
   - Jon reviews and edits every sentence and checks every reference.
   - Validation (null and injection calibration, gaiasupdate reproduction, Gaia-4b recovery) is documented in `notes/phase1_validation.md`.
3. **References:** export BibTeX from ADS, or verify each entry against ADS/arXiv by hand. Never unverified entries.
4. **Figure:** if the Note uses a figure, its caption names the tool if Claude produced the plotting code. The one-table-or-one-figure limit applies.
5. **The duotransits note (submitted 2026-10-12):** check whether it carries an AI disclosure statement. If the editor asks, add one in the same form.

## RNAAS preparation rules

Source: https://journals.aas.org/research-note-preparation-guidelines

**Length**
- At most 1,500 words, counting:
  - title, headers, captions and references;
  - acknowledgments of people, facilities and software;
  - about 150 words reserved for the required abstract.
- Table contents do not count; captions do.
- The editorial office's word count is final.

**Content**
- One figure **or** one table, not both. A figure plus its caption must fit one page.
- Long tables show 5–10 lines in print; the full table is machine-readable (MRT) or on Zenodo (AAS community).
- The data behind a figure are encouraged.

**Format**
- No line numbers. AASTeX v7+ preferred; Overleaf accepted.
- Single author writes in the first person singular.

**Process**
- Moderated, not peer reviewed; typically published within 72 h.
- Published as submitted, with no copy-editing.
- DOI and ADS indexing; CC-BY 4.0; no charges.
- At most 20 Notes per first author per 12 months.

## Lessons from the submitted duotransits note (AAS82207, 2026-10-12)

Reported by the duo-transits session. Submitted source: `C:\git\duo-transits`, `singles/rnaas/rnaas.tex` at commit d932e46.

**Word limit:** the submission site states **1,350 words** plus a 150-word abstract, with one figure **or** one table. The guidelines page says 1,500 including references. Plan for 1,350.

**Template**
- Overleaf AAS template "AASTeX Template for submissions to AAS Journals (ApJ-AJ-ApJS-ApJL-PSJ-RNAAS)".
- Class line: plain `\documentclass{aastex701}`. The RNAAS option is deprecated in v7.
- Author block: `\author[orcid=0009-0000-7524-4470]{Jon McBee}`, `\affiliation{Independent researcher, USA}`, `\email{...}`. Jon confirmed this form on 2026-10-12.

**AI disclosure**
- A `\section*{AI disclosure statement}` after `\facilities`/`\software` and before the bibliography.
- The submission site has no separate AI field.

**References**
- That note used a hand-written `thebibliography`, with every entry verified against arXiv/ADS/ASCL metadata.
- Either that or ADS-exported BibTeX meets the guideline. Never unverified entries.

**Submission route**
- Overleaf "Submit to AAS" button, then an emailed link to aas.msubmit.net.
- The site sits behind a Cloudflare check, so Jon must finish this step in his own browser.

**On the site**
- **Corridor:** "The Solar System, Exoplanets, and Astrobiology".
- **orcid-ID.png:** retype it from "Figure file" to "Other LaTeX file", or the note appears to have a figure as well as a table.
- **Abstract field:** the imported text contains raw LaTeX (`\,`, `$...$`, `\%`). Replace it with plain text.
- **Dataset DOIs:** answering Yes for MAST requires DOIs; Yes for IPAC asks for citations instead. For this note: NASA Exoplanet Archive PSCompPars 10.26133/NEA13 (verified on DataCite) and the Gaia DR4 DOI.
- **Other answers:** digital assets none; funding no; ethics statement yes.
- **UAT concepts:** entered on their own page.
- **Approval:** the Approve checkbox appears only after you open the merged PDF.
