"""Independent reconstruction of the frozen TN10 four-cell contrast.

Does not import card002_run. Uses a different PDF extractor and cell-regression
HC2 covariance as a cross-check of the direct variance formula. This is a
second implementation on the same data.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--raw-root', type=Path, help='Read-only original input directory')
    args = parser.parse_args()
    root = args.data_root.resolve()
    raw_root = args.raw_root.resolve() if args.raw_root is not None else root / 'raw'
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    from scipy.stats import t
    from pypdf import PdfReader

    table = []
    for page in PdfReader(raw_root / 'JCI176403_Table6.pdf').pages:
        for line in (page.extract_text() or '').splitlines():
            fields = line.split()
            if len(fields) != 3:
                continue
            try:
                correlation, pvalue = float(fields[1]), float(fields[2])
            except ValueError:
                continue
            if correlation < 0:
                table.append(fields[0])
    assert len(table) == len(set(table)) == 359
    meta = pd.read_csv(raw_root / 'TN10_metadata.csv')
    counts = pd.read_csv(raw_root / 'TN10_raw_counts.csv', index_col=0)
    base = meta.query('`Time.Point` == "Baseline"').copy()
    symbols = sorted(set(table).intersection(counts.index))
    digest = hashlib.sha256(('\n'.join(symbols)+'\n').encode()).hexdigest()
    assert digest == '09c7b17296f6ec5b1b90357b555ec62748821592f042070a6f0ec525cd106f11'
    # Data-frame axis operations deliberately differ from production arrays.
    raw = counts.loc[:, base.libraryID]
    logcpm = np.log2(raw.loc[symbols].div(raw.sum(axis=0), axis=1) * 1e6 + 1)
    standardized = logcpm.sub(logcpm.mean(axis=1), axis=0).div(logcpm.std(axis=1, ddof=1), axis=0)
    score = standardized.mean(axis=0)
    high = (score > score.median()).to_numpy()
    arm = base.TREATMENTDESC.eq('Teplizumab').to_numpy()
    outcome = np.minimum(base.Survival.to_numpy(float), 156)
    cell_masks = [arm & ~high, ~arm & ~high, arm & high, ~arm & high]
    cells = [outcome[mask] for mask in cell_masks]
    signs = np.array([1, -1, -1, 1])
    n = np.array([len(c) for c in cells])
    means = np.array([c.mean() for c in cells])
    terms = np.array([c.var(ddof=1)/len(c) for c in cells])
    theta = float(signs @ means)
    se = math.sqrt(terms.sum())
    df = float(terms.sum()**2 / np.sum(terms**2/(n-1)))
    ci = (theta + np.array([-1, 1])*t.ppf(.975, df)*se).tolist()
    design = np.column_stack(cell_masks).astype(float)
    model = sm.OLS(outcome, design).fit(cov_type='HC2')
    regression_theta = float(signs @ model.params)
    regression_se = float(np.sqrt(signs @ model.cov_params() @ signs))
    expected = json.loads((root / 'results' / 'primary_results.json').read_text())
    checks = {
        'independent_pdf_gene_projection': digest,
        'cell_counts': n.tolist(),
        'cell_means_weeks': means.tolist(),
        'theta_weeks': theta,
        'se_weeks': se,
        'welch_df': df,
        'ci95_weeks': ci,
        'hc2_regression_theta': regression_theta,
        'hc2_regression_se': regression_se,
        'matches_primary': bool(np.allclose([theta, se, *ci], [expected['theta']['estimate_weeks'], expected['theta']['se_weeks'], *expected['theta']['ci95_weeks']], rtol=0, atol=1e-10)),
        'matches_independent_regression': bool(np.allclose([theta, se], [regression_theta, regression_se], rtol=0, atol=1e-10)),
        'scope': 'Independent code reconstruction on identical public data; not external validation',
    }
    assert checks['matches_primary'] and checks['matches_independent_regression']
    (root / 'results' / 'independent_verification.json').write_text(json.dumps(checks, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
