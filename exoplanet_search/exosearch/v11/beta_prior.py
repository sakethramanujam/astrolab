"""Write configs/v1.1_beta_prior.json: the empirical-Bayes prior of beta option d (beta_treatment.md), per sector and
trial duration: median and between-star variance of the star-quality beta_D of kept stars (exosearch.stat_options.prior,
star_quality_v1 tables; no null, injection or candidate row enters). Shipped as a file so every host uses the same
numbers; its sha256 is part of the v1.1 config hash.
  python -m exosearch.v11.beta_prior
"""
from __future__ import annotations

import json

from .. import quality, stat_options
from . import settings


def main(sectors=(70, 106)):
    q = quality.load(sectors)
    out = {str(s): {k: [float(v[0]), float(v[1])] for k, v in stat_options.prior(q[q.sector == s]).items()} for s in sectors}
    out["_note"] = "per sector: {duration: [median beta_D, tau^2]} from star_quality_v1 kept stars; exosearch.v11.beta_prior"
    p = settings.ROOT / settings.cfg()["beta_prior_file"]
    p.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(p.read_text())


if __name__ == "__main__":
    main()
