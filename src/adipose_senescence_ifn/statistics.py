"""Animal-label permutation inference and multiplicity control."""

from __future__ import annotations

from math import factorial

import numpy as np
import pandas as pd


def _slope(age: np.ndarray, outcome: np.ndarray) -> float:
    design = np.column_stack([np.ones(len(age)), age])
    coefficient, *_ = np.linalg.lstsq(design, outcome, rcond=None)
    return float(coefficient[1])


def _unique_multiset_permutations(values: list[float]):
    counts: dict[float, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1

    def generate(prefix: list[float]):
        if len(prefix) == len(values):
            yield tuple(prefix)
            return
        for value in sorted(counts):
            if counts[value] == 0:
                continue
            counts[value] -= 1
            prefix.append(value)
            yield from generate(prefix)
            prefix.pop()
            counts[value] += 1

    yield from generate([])


def permutation_slope(
    frame: pd.DataFrame,
    outcome: str,
    n_permutations: int = 9999,
    seed: int = 724992,
) -> dict:
    required = {"animal_id", "age_months", outcome}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing columns: {sorted(required - set(frame.columns))}")
    animal_age = frame[["animal_id", "age_months"]].drop_duplicates()
    if animal_age["animal_id"].duplicated().any():
        raise ValueError("Age must be constant within animal")
    if animal_age["animal_id"].nunique() < 3 or animal_age["age_months"].nunique() < 2:
        raise ValueError("At least three animals and two ages are required")

    animals = animal_age["animal_id"].tolist()
    ages = animal_age["age_months"].astype(float).tolist()
    observed = _slope(
        frame["age_months"].to_numpy(dtype=float),
        frame[outcome].to_numpy(dtype=float),
    )

    denominator = 1
    for count in pd.Series(ages).value_counts():
        denominator *= factorial(int(count))
    unique_count = factorial(len(ages)) // denominator

    if unique_count <= n_permutations + 1:
        permutations = _unique_multiset_permutations(ages)
        method = f"exact ({unique_count} unique assignments)"
    else:
        rng = np.random.default_rng(seed)
        permutations = (tuple(rng.permutation(ages)) for _ in range(n_permutations))
        method = f"Monte Carlo ({n_permutations} assignments; seed={seed})"

    lookup_rows = frame["animal_id"].map({animal: i for i, animal in enumerate(animals)})
    permuted_slopes: list[float] = []
    for permuted in permutations:
        permuted_age = np.asarray(permuted, dtype=float)[lookup_rows.to_numpy()]
        permuted_slopes.append(
            _slope(permuted_age, frame[outcome].to_numpy(dtype=float))
        )
    exceedances = sum(abs(value) >= abs(observed) for value in permuted_slopes)
    if unique_count <= n_permutations + 1:
        p_value = exceedances / len(permuted_slopes)
    else:
        # The +1 correction prevents zero Monte Carlo P values.
        p_value = (exceedances + 1) / (len(permuted_slopes) + 1)

    youngest = float(frame["age_months"].min())
    oldest = float(frame["age_months"].max())
    contrast = float(
        frame.loc[frame["age_months"] == oldest, outcome].mean()
        - frame.loc[frame["age_months"] == youngest, outcome].mean()
    )
    return {
        "slope_per_month": observed,
        "oldest_minus_youngest": contrast,
        "permutation_p": p_value,
        "permutation_method": method,
        "n_animals": int(frame["animal_id"].nunique()),
        "n_rows": len(frame),
    }


def bootstrap_slope_interval(
    frame: pd.DataFrame,
    outcome: str,
    n_bootstrap: int = 2000,
    seed: int = 724992,
) -> tuple[float, float]:
    animals = frame["animal_id"].unique()
    rng = np.random.default_rng(seed)
    slopes: list[float] = []
    for _ in range(n_bootstrap):
        sampled = rng.choice(animals, size=len(animals), replace=True)
        pieces = []
        for draw, animal in enumerate(sampled):
            piece = frame[frame["animal_id"] == animal].copy()
            piece["animal_id"] = f"draw-{draw}"
            pieces.append(piece)
        resample = pd.concat(pieces, ignore_index=True)
        slopes.append(
            _slope(
                resample["age_months"].to_numpy(dtype=float),
                resample[outcome].to_numpy(dtype=float),
            )
        )
    low, high = np.quantile(slopes, [0.025, 0.975])
    return float(low), float(high)


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    if values.empty:
        return values.copy()
    order = np.argsort(values.to_numpy(dtype=float))
    ranked = values.to_numpy(dtype=float)[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=values.index)


def run_inference(
    aggregates: pd.DataFrame,
    score_columns: list[str],
    n_permutations: int = 9999,
    n_bootstrap: int = 2000,
) -> pd.DataFrame:
    results: list[dict] = []
    for (depot, subtype), group in aggregates.groupby(["depot", "subtype"]):
        if group["animal_id"].nunique() < 3 or group["age_months"].nunique() < 2:
            continue
        for score in score_columns:
            estimate = permutation_slope(
                group, score, n_permutations=n_permutations
            )
            low, high = bootstrap_slope_interval(
                group, score, n_bootstrap=n_bootstrap
            )
            results.append(
                {
                    "depot": depot,
                    "subtype": subtype,
                    "score": score,
                    **estimate,
                    "slope_ci_low": low,
                    "slope_ci_high": high,
                }
            )
    result = pd.DataFrame(results)
    if not result.empty:
        result["q_value"] = (
            result.groupby("score", group_keys=False)["permutation_p"]
            .apply(benjamini_hochberg)
            .sort_index()
        )
    return result
