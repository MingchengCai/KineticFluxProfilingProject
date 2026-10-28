# SP-KFP flux calculator

A standardized, GUI-based data-processing pipeline for **non-stationary isotope
tracing / kinetic flux profiling (SP-KFP)**. It reads the .xlsx labeling
template, fits one `k_app` per metabolite, lets you choose the pseudo tracer and
the metabolites to pool, checks Criterion 1 and reports the flux with its
standard error.

```
python spkfp_app.py                                     # start the GUI
python spkfp_app.py --check data.xlsx                   # audit + k_app table
python spkfp_app.py --check data.xlsx --pt IMP --flux AMP,ADP,ATP
```

---

## 1. Requirements — installed automatically

| Component | Needed | Notes |
|---|---|---|
| Python | 3.8 – 3.12 | the only thing that must already be on the machine |
| tkinter | yes | part of Python, but an **optional component of the installer** |
| numpy | ≥ 1.20 | **installed on first run if missing** |
| matplotlib | ≥ 3.4 | **installed on first run if missing** |
| openpyxl | ≥ 3.0 | **installed on first run if missing** |

### What has to be downloaded, and what does not

Everything the program imports falls into exactly three groups:

| Group | Modules | Where they come from | In the auto-install? |
|---|---|---|---|
| **Standard library** | `os`, `sys`, `math`, `re`, `subprocess`, `tempfile`, `traceback`, `importlib.util`, `dataclasses`, `typing` | shipped inside every CPython installation | **no** — nothing to download |
| **Standard library, optional component** | `tkinter` (`tkinter.ttk`, `tkinter.filedialog`, `tkinter.messagebox`) | ships with the python.org / Anaconda installers, but is a **tick box** ("tcl/tk and IDLE") on Windows and a **separate distro package** on Linux (`python3-tk`) | **cannot** — pip has no `tkinter` |
| **Third-party** | `numpy`, `matplotlib`, `openpyxl` | site-packages; must be downloaded once | **yes**, exactly these three |

So on a freshly installed Python nothing needs to be fetched except the three
third-party packages, which the program fetches itself. (`pip install os` is not
even a thing — the standard library is not on PyPI.)

The only thing the program cannot repair is a Python that was installed
**without tcl/tk**. It then shows a native message box (no tkinter needed to
display it) explaining what to do:

* **Windows** — run the Python installer again → *Modify* → tick
  *tcl/tk and IDLE*;
* **macOS** — use the installer from python.org (it includes tkinter);
* **Linux** — `sudo apt install python3-tk` (or the distribution equivalent).

`--check` (command line, no window) works even without tkinter.

**You can hand this file to a colleague whose computer only has Python.** At
start-up the three third-party packages are imported in a guarded way; if one is
missing, the program

1. shows a "please wait" window,
2. runs `"<this python>" -m pip install numpy matplotlib openpyxl` (trying
   PyPI first, then the Tsinghua and Aliyun mirrors, and `--user` if the
   environment is not writable; `pip` itself is enabled with `ensurepip` when
   absent),
3. restarts itself automatically, now fully equipped.

If the automatic installation is impossible (no internet, proxy, locked-down
machine) an error box shows the *exact* command to run once by hand, with the
full path of the interpreter.

### Using the program outside China

The package indexes are tried in this order:

1. **your own index**, if `SPKFP_PIP_INDEX_URL` is set;
2. **PyPI** — the default of pip, so a `pip.conf` / `%APPDATA%\pip\pip.ini`
   pointing at a company mirror is honoured here;
3. the **Tsinghua** and **Aliyun** mirrors, only as a fallback.

So on a computer outside China the normal path is step 2 and the Chinese
mirrors are never contacted.  They are public HTTPS `simple` indexes (no
account, no geo-restriction) and do work from abroad, only with higher latency,
so keeping them as a last resort does no harm.

Two safeguards keep a bad network from blocking the program:

* every index is **probed for 5 s** before pip is started on it, so a dead or
  firewalled index is skipped immediately instead of after a long pip timeout
  (an offline computer therefore fails in seconds with a clear message);
* an attempt that starts but is **too slow is abandoned** after 7 min and the
  next index is tried, with a hint to set `SPKFP_PIP_INDEX_URL` if a different
  mirror would be faster.

If you prefer to drop the Chinese mirrors entirely, or to fix the order
yourself:

| Variable | Effect |
|---|---|
| `SPKFP_NO_AUTO_INSTALL=1` | never install anything; report what is missing instead |
| `SPKFP_SKIP_CN_MIRRORS=1` | never fall back to the Chinese mirrors |
| `SPKFP_PIP_INDEX_URL=<url>` | try this index first (e.g. a company mirror) |
| `SPKFP_PIP_INDEXES=<url1>;<url2>` | use exactly this chain, in this order |
| `SPKFP_PIP_TIMEOUT=<seconds>` | abandon one index after this long (default 300) |
| `HTTPS_PROXY` / `HTTP_PROXY` | pip honours them, so corporate proxies work as usual |
| `SPKFP_BOOTSTRAP_DONE` | set internally after a restart, do not set it yourself |

No scipy / pandas / R is required — the fits are implemented in numpy.  (scipy
is used by the test suite only, for an independent cross-check.)


---

## 2. Workflow

| Step | What happens |
|---|---|
| **1. Submit data** | `Load data (.xlsx)` opens a file dialog. The workbook is audited; if something is wrong a dialog lists **every** problem and the file dialog re-opens so a corrected file can be submitted. |
| **2. Fit every metabolite** | Each labeled form is normalized by the mean labeling fraction of its **last time point**, then `L(t) = 1 − exp(−k_app·t)` is fitted with `k_app` as the only parameter. The table shows `k_app ± SE`, R², n and the normalization reference; the figure shows all normalized curves. `Save figure as image …` writes a 300 dpi PNG/TIFF/PDF/SVG. |
| **3. Choose the flux inputs** | `Calculate flux` opens a dialog: **one** pseudo tracer (combobox) and **one or more** metabolites to pool (checkboxes). The dialog states that the pooled metabolites must be **downstream** of the pseudo tracer, i.e. that Criterion 3 (major branch) and Criterion 2 (reversibility) hold for this group. Metabolites without a pool size are shown but cannot be ticked. |
| **4. Criterion 1** | `k_app(pseudo tracer) ≥ 4 × k_app` of **every** pooled metabolite. On failure a dialog names each offending pair with both values and the threshold; the flux is not calculated. |
| **5. Flux** | The labeled concentrations of the pooled metabolites are summed per sample and fitted with `L = Plateau·(1 − exp(−k_app·t))`; `flux = k_app × Plateau ± SE`. |
| **6. Output** | The flux panel shows the pseudo tracer, the pooled metabolites, `k_app`, `Plateau`, R², n and the Criterion-1 comparison lines. `Export results …` writes an .xlsx with the k_app table, the Criterion-1 report and the flux. |

---

## 3. Input workbook

### Sheet 1 — labeling fractions (first worksheet)

Minimal layout (the column of metabolite names may be **without a header**):

| Metabolite | isotopeLabel | formula | 0_1 | 0_2 | … | 6_3 |
|---|---|---|---|---|---|---|
| Hypoxanthine | C12 PARENT | C5H4N4O | 1 | 1 | | 0.05 |
| Hypoxanthine | N15-label-4 | C5H4N4O | 0 | 0 | | 0.95 |

* **Header row** — the row must contain `isotopeLabel` and `formula`.
  `Metabolite` and `Group` headers are optional:
  * with a `Metabolite` header that column is used;
  * without it, the names are read from the column **immediately left of**
    `isotopeLabel` (as in the supplied template);
  * a `Group` column, if present, is ignored (informational only).
* **Sample columns** — everything to the right of `formula`, named
  `<time>_<sample>` (e.g. `0_1`, `0.05_2`, `6_3`). The part before `_` is the
  time after tracer addition in hours, the part after `_` is the replicate.
  At least two different time points are required; the number of replicates may
  differ between time points.
* **Rows** — one per metabolite × isotopologue. Every row needs the metabolite
  name and the isotope label. `PARENT` marks the unlabeled isotopologue;
  **every metabolite needs at least one non-PARENT row**. Values are
  proportions in 0 … 1 (percentage data are rejected with a hint).
* If a metabolite has more than one labeled form, they are averaged with
  `SUM(fraction × labeled atoms) / (atoms of that element)` read from the label
  and the formula. With a single labeled form the raw fraction is used.

### Sheet 2 — pool sizes (second worksheet)

```
A1: (empty)   B1: pmole/million cells
A2: AMP       B2: 107.507
A3: ADP       B3: 850.089
A4: ATP       B4: 13141.756
```

* `B1` is the unit (the flux unit is that unit per hour).
* Column A = metabolite name, column B = pool size (a number; `0` allowed).
* A pool size is **only needed for the metabolites selected in the flux
  dialog**. Metabolites without one can still be fitted, they simply cannot be
  pooled; the load note lists them.

---

## 4. Audits performed on submission

The whole list of problems is reported in one dialog:

* two worksheets exist;
* the header row contains `isotopeLabel` and `formula`;
* every column right of `formula` is a valid `<time>_<sample>` header, no
  duplicates, at least two time points;
* every data row carries a metabolite name and an isotope label;
* every metabolite has at least one non-PARENT row;
* values are numeric and inside 0–1;
* pool sizes (where given) are numeric and not negative.

Blank cells are tolerated (they are dropped from the fit and reported as a
note).

---

## 5. Model and estimators

### Per-metabolite k_app

The labeling fraction of the labeled form is divided by its mean at the **last
time point** (`L_last`), so that `L = 1` there, and

```
L(t) = 1 − exp(−k_app · t)
```

is fitted by least squares on all individual samples (replicates independent).
The only parameter is `k_app`; there is no `Plateau`.

```
SE(k_app) = sqrt( RSS / (n−1) / Σ (t·exp(−k_app·t))² )
```

If the labeling is already complete at the first sampled time point
(`k_app·t_min ≥ 5`) the curve is flat at 1 and `k_app` is **not identifiable**;
the table then shows `n.i. (inst.)` and Criterion 1 treats such a metabolite as
infinitely fast (a pseudo tracer like that always passes; using one as a pooled
metabolite always fails, since nothing can be 4× faster).

The fit is solved by scanning the residual sum of squares over a wide log-`k`
range and refining with a golden-section search — deterministic, no starting
values, agrees with `scipy.optimize.curve_fit` to ~1e-6.

### Flux

```
L_pooled(t) = Σ labeled concentration of the selected metabolites
L_pooled(t) = Plateau · (1 − exp(−k_app · t))

flux    = k_app · Plateau
SE²     = Plateau²·SE(k_app)² + k_app²·SE(Plateau)²
```

`k_app` and `Plateau` are treated as **independent** parameters, so their
covariance is **not** propagated.  (The least-squares estimates of the two
parameters are in fact strongly correlated — on the template data
ρ ≈ −0.996 — so this choice gives a deliberately conservative, larger standard
error; the fitted covariance is still computed and kept in
`FitResult.cov_kp` / `FluxResult.cov_kp`, it is just not used in the SE.)

### Criterion 1

```
k_app(pseudo tracer) ≥ 4 × k_app(every pooled metabolite)
```

is checked for each pooled metabolite after the dialog, and the comparison
lines are shown in the flux panel and saved in the export.

---

## 6. Files

```
spkfp_app.py              the application (GUI + analysis + --check CLI)
tests/test_spkfp.py       automated tests (55 assertions) on modified templates
screenshots/              gui_loaded.png, gui_flux.png, dialog_select.png,
                          fitted_curves.png
README.md                 this document
```

![GUI after loading a file](screenshots/gui_loaded.png)

![Flux selection dialog](screenshots/dialog_select.png)

![Saved figure](screenshots/fitted_curves.png)

Run the tests (they build modified copies of the template in a temp folder):

```
D:\Anaconda3\python.exe tests\test_spkfp.py "D:\Desktop\Template for SP-KFP.xlsx"
```

Expected output: `56 passed, 0 failed`.

---

## 7. Reference result for the supplied template

```
metabolite      k_app ± SE (h⁻¹)      R²      L_last
Hypoxanthine    n.i. (instantaneous) 0.995   0.9484
Inosine         6.088 ± 1.195        0.789   0.8824
IMP             2.151 ± 0.366        0.843   0.9251
AMP             0.2613 ± 0.0329      0.897   0.5642
ADP             0.3121 ± 0.0273      0.950   0.6034
ATP             0.3299 ± 0.0264      0.959   0.6014
```

With pseudo tracer **IMP** and pooled metabolites **AMP, ADP, ATP**:

```
Criterion 1: PASS
  IMP 2.151 h⁻¹ ≥ 4 × 0.2613 (AMP, ratio 8.23)
  IMP 2.151 h⁻¹ ≥ 4 × 0.3121 (ADP, ratio 6.89)
  IMP 2.151 h⁻¹ ≥ 4 × 0.3299 (ATP, ratio 6.52)

k_app (pooled) = 0.1032 ± 0.0206 h⁻¹
Plateau        = 1.855e4 ± 2.84e3 pmole/10⁶ cells
Flux           = 1915 ± 481.5 pmole/10⁶ cells/h   (R² = 0.994, n = 21)
```

Note: the fitted Plateau equals 132 % of the pooled pool size, i.e. the labeling
had not reached its plateau inside the 6 h window and the flux relies on
extrapolation — the panel says so explicitly.
