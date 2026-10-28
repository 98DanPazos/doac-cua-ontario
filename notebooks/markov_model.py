"""
markov_model.py
Cohort Markov model: generic apixaban vs warfarin in atrial fibrillation.

5 states, 3-month cycles, age 70 to 100, CADTH worked-example structure.
Every parameter can be a single number (deterministic run) or a NumPy array
of draws (probabilistic run, Step 6). The same code handles both.
"""
import json
import numpy as np
import pandas as pd

STATES = ["S0", "S1", "S2", "S3", "Dead"]
REF_AGE = 70          # age at which the trial-based rates apply


# --- conversion helpers (Step 3, part 1) -----------------------------
def to_rate(prob):
    return -np.log(1 - prob)

def to_cycle_prob(rate, dt):
    return 1 - np.exp(-rate * dt)

def apply_odds_ratio(prob, odds_ratio):
    return odds_ratio * prob / (1 - prob + odds_ratio * prob)


def load_inputs(folder="../outputs"):
    """Return (p, settings): base-case parameters and model settings."""
    inputs = pd.read_csv(f"{folder}/inputs.csv", index_col="name")
    with open(f"{folder}/settings.json") as f:
        settings = json.load(f)
    return inputs["base"].to_dict(), settings


def run_arm(p, arm, settings, keep_trace=False):
    """
    Run one treatment arm through the full lifetime horizon.

    p         dict of parameters (numbers or arrays)
    arm       "apixaban" or "warfarin"
    settings  dict from settings.json (can be edited for scenarios)
    Returns a dict of discounted and undiscounted totals.
    """
    dt = settings["cycle_years"]
    n_cycles = int(round((settings["max_age"] - settings["start_age"]) / dt))
    slope, intercept = settings["gompertz_slope"], settings["gompertz_intercept"]

    def qx_general(age):
        # Fitted StatCan mortality at a given age (Gompertz curve from Step 2)
        return np.minimum(np.exp(intercept + slope * age), 1.0)

    # Everyone starts in S0. Each state holds the FRACTION of the cohort in it.
    occ = {s: 0.0 for s in STATES}
    occ["S0"] = 1.0

    # Annual cost of the strategy while in S0: the ONLY drug difference
    if arm == "apixaban":
        cost_S0_drug = p["c_apixaban"]
    else:
        cost_S0_drug = p["c_warfarin"] + p["c_inr"]

    # Annual utility and long-term cost of each living state
    u_major = (1 - p["p_rankin_5"]) * p["u_rankin_3_4"] + p["p_rankin_5"] * p["u_rankin_5"]
    utility = {"S0": p["u_af"], "S1": p["u_minor_stroke"], "S2": u_major, "S3": p["u_af"]}
    lt_cost = {"S0": 0.0, "S1": p["c_lt_minor_stroke"], "S2": p["c_lt_major_stroke"], "S3": 0.0}

    # Running totals
    keys = ["cost", "cost_drug", "cost_event", "cost_lt", "qaly", "ly",
            "cost_undisc", "qaly_undisc", "ly_undisc", "strokes", "ich", "major_bleeds"]
    tot = {k: 0.0 for k in keys}
    trace, yearly_cost, year_acc = [], [], 0.0

    for t in range(n_cycles):
        age = settings["start_age"] + t * dt
        discount = (1 + settings["discount_rate"]) ** -(t * dt)

        # Optional scenario: apixaban's benefit fades after `wane_start` years.
        # eff = 1 means full effect; eff -> 0 means OR -> 1 (no effect).
        # Defaults (no waning) apply when the keys aren't in settings.
        eff = 1.0
        if t * dt >= settings.get("wane_start", 999):
            eff = (1 - settings.get("wane_decay", 0)) ** (t * dt - settings["wane_start"])

        # Background death: the trial's rate at 70, rising with age like
        # the general population (only the GRADIENT of the life table is used)
        q_bg = np.minimum(p["p_nonvasc_death"] * qx_general(np.floor(age))
                          / qx_general(REF_AGE), 0.999)

        if keep_trace:
            trace.append({"age": age, **{s: float(np.mean(occ[s])) for s in STATES}})

        new = {s: 0.0 for s in STATES}
        new["Dead"] = occ["Dead"]
        cycle_cost = 0.0

        for s in ["S0", "S1", "S2", "S3"]:
            prior_stroke = s in ("S1", "S2")        # stroke history raises risks
            on_aspirin = s != "S0"                  # anticoagulation stopped

            # 1. background death (higher after a stroke)
            hr = p["hr_death_prior_stroke"] if prior_stroke else 1.0
            p_death = to_cycle_prob(to_rate(q_bg) * hr, dt)
            alive = occ[s] * (1 - p_death)

            # 2. event RATES, adjusted for age, history and aspirin
            r = {
                "stroke": to_rate(p["p_stroke"]) * p["rr_stroke_per_decade"] ** ((age - REF_AGE) / 10),
                "ich":    to_rate(p["p_ich"]),
                "mb":     to_rate(p["p_major_bleed"]),
                "mi":     to_rate(p["p_mi"]),
                "tia":    to_rate(p["p_tia"]),
                "minor":  to_rate(p["p_minor_bleed"]),
            }
            if prior_stroke:
                r["stroke"] = r["stroke"] * p["rr_stroke_prior"]
            if on_aspirin:
                for ev, rr in [("stroke", "rr_asa_stroke"), ("ich", "rr_asa_ich"),
                               ("mb", "rr_asa_mb"), ("mi", "rr_asa_mi"),
                               ("tia", "rr_asa_tia"), ("minor", "rr_asa_minor")]:
                    r[ev] = r[ev] * p[rr]

            prob = {ev: to_cycle_prob(rate, dt) for ev, rate in r.items()}
            if arm == "apixaban" and s == "S0":
                for ev, orn in [("stroke", "or_apx_stroke"), ("ich", "or_apx_ich"),
                                ("mb", "or_apx_mb"), ("mi", "or_apx_mi")]:
                    prob[ev] = apply_odds_ratio(prob[ev], p[orn] ** eff)

            n = {ev: alive * pr for ev, pr in prob.items()}

            # 3. route events
            f_stroke = np.minimum(p["f_stroke_fatal"] * (p["hr_recur_fatal"] if prior_stroke else 1.0), 1.0)
            stroke_dead  = n["stroke"] * f_stroke
            stroke_major = n["stroke"] * (1 - f_stroke) * p["f_stroke_major"]
            stroke_minor = n["stroke"] * (1 - f_stroke) * (1 - p["f_stroke_major"])
            ich_live = n["ich"] * (1 - p["f_bleed_fatal"])
            mb_live  = n["mb"] * (1 - p["f_bleed_fatal"])
            mi_dead  = n["mi"] * p["f_mi_fatal"]
            mi_live  = n["mi"] - mi_dead
            stay = alive - n["stroke"] - n["ich"] - n["mb"] - mi_dead

            new["Dead"] += (occ[s] * p_death + stroke_dead + (n["ich"] - ich_live)
                            + (n["mb"] - mb_live) + mi_dead)
            new["S2"] += stroke_major
            # minor stroke / ICH: to S1, unless already in S2 (stay there)
            new["S2" if s == "S2" else "S1"] += stroke_minor + ich_live
            # major bleed: S0 -> S3; already on aspirin -> stay put
            new["S3" if s == "S0" else s] += mb_live
            new[s] += stay

            # 4. costs: drug + long-term care on time in state, plus event costs
            c_drug = occ[s] * (cost_S0_drug if s == "S0" else p["c_aspirin"]) * dt
            c_lt = occ[s] * lt_cost[s] * dt
            c_evt = (stroke_dead * p["c_stroke_fatal"] + stroke_major * p["c_stroke_major"]
                     + stroke_minor * p["c_stroke_minor"] + n["ich"] * p["c_ich"]
                     + n["mb"] * p["c_major_bleed"] + mi_dead * p["c_mi_fatal"]
                     + mi_live * p["c_mi_nonfatal"] + n["tia"] * p["c_tia"]
                     + n["minor"] * p["c_minor_bleed"])

            # 5. QALYs: time in state x utility, minus one-off event losses
            u = utility[s] - p["dec_per_year_over_70"] * max(age - REF_AGE, 0)
            q = occ[s] * u * dt - dt * (mi_live * p["dec_mi"] + mb_live * p["dec_major_bleed"]
                                        + n["minor"] * p["dec_minor_bleed"] + n["tia"] * p["dec_tia"])

            # 6. add to totals (discounted and undiscounted)
            tot["cost"] += (c_drug + c_lt + c_evt) * discount
            tot["cost_drug"] += c_drug * discount
            tot["cost_event"] += c_evt * discount
            tot["cost_lt"] += c_lt * discount
            tot["qaly"] += q * discount
            tot["ly"] += occ[s] * dt * discount
            tot["cost_undisc"] += c_drug + c_lt + c_evt
            tot["qaly_undisc"] += q
            tot["ly_undisc"] += occ[s] * dt
            tot["strokes"] += n["stroke"]
            tot["ich"] += n["ich"]
            tot["major_bleeds"] += n["mb"]
            cycle_cost = cycle_cost + c_drug + c_lt + c_evt

        occ = new
        year_acc = year_acc + cycle_cost
        if (t + 1) % 4 == 0:                    # every 4 cycles = 1 year
            yearly_cost.append(year_acc)
            year_acc = 0.0

    tot["yearly_cost_undisc"] = yearly_cost     # used by the budget impact (Step 8)
    if keep_trace:
        tot["trace"] = pd.DataFrame(trace)
    return tot


def run_both(p, settings, **kw):
    """Run both arms with the same inputs. Returns (apixaban, warfarin)."""
    return run_arm(p, "apixaban", settings, **kw), run_arm(p, "warfarin", settings, **kw)
