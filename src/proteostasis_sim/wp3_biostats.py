from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import minimize


class CompetingRisksBiostats:
    """WP3: Shared-frailty competing risks, exact Aalen-Johansen, Cox CSH, Fine-Gray IPCW, and Phi index."""

    def __init__(self, interp_gen, interp_amb, shape_gen: float = 3.2, shape_amb: float = 2.5):
        self.interp_gen = interp_gen
        self.interp_amb = interp_amb
        self.shape_gen = float(shape_gen)
        self.shape_amb = float(shape_amb)

    def simulate_virtual_cohort(
        self,
        E_vals: np.ndarray,
        n_per_arm: int = 2000,
        mode: str = "hormese",
        sigma_shared: float = 0.25,
        sigma_ind: float = 0.22,
        t_cens: float = 150.0,
        random_censor_rate: float = 0.0,
        seed: int = 42,
    ) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        sc_g0 = float(self.interp_gen(0.0))
        blocks = []

        for E_val in E_vals:
            sc_a = float(np.clip(self.interp_amb(E_val), 5.0, 2500.0))
            sc_g = float(np.clip(self.interp_gen(E_val), 5.0, 2500.0)) if mode == "hormese" else sc_g0

            z_shared = rng.normal(0.0, sigma_shared, size=n_per_arm)
            norm_shift = -0.5 * (sigma_shared**2 + sigma_ind**2)
            theta_g = np.exp(z_shared + rng.normal(0.0, sigma_ind, size=n_per_arm) + norm_shift)
            theta_a = np.exp(z_shared + rng.normal(0.0, sigma_ind, size=n_per_arm) + norm_shift)

            T_gen = sc_g * (rng.weibull(self.shape_gen, size=n_per_arm) / theta_g)
            T_amb = sc_a * (rng.weibull(self.shape_amb, size=n_per_arm) / theta_a)
            T_chr = 0.55 * T_amb

            if random_censor_rate > 0.0:
                T_c = np.minimum(rng.exponential(1.0 / random_censor_rate, size=n_per_arm), t_cens)
            else:
                T_c = np.full(n_per_arm, t_cens, dtype=float)

            T_min = np.minimum(T_gen, T_amb)
            T_obs = np.minimum(T_min, T_c)
            event = np.zeros(n_per_arm, dtype=int)
            occ = T_min <= T_c
            event[occ & (T_gen <= T_amb)] = 1
            event[occ & (T_amb < T_gen)] = 2

            blocks.append(pd.DataFrame({
                "E": float(E_val),
                "time": T_obs,
                "event": event,
                "T_gen_lat": T_gen,
                "T_amb_lat": T_amb,
                "T_chr_lat": T_chr,
            }))
        return pd.concat(blocks, ignore_index=True)

    @staticmethod
    def aalen_johansen_cif(times: np.ndarray, events: np.ndarray, target_cause: int, t_ref: float = 100.0) -> float:
        times = np.asarray(times, dtype=float)
        events = np.asarray(events, dtype=int)
        order = np.argsort(times)
        t_s, e_s = times[order], events[order]
        unique_t = np.unique(t_s[e_s > 0])
        S_prev = 1.0
        t_pts, c_pts = [0.0], [0.0]
        for ut in unique_t:
            at_risk = int(np.sum(t_s >= ut))
            if at_risk == 0:
                break
            d_all = int(np.sum((t_s == ut) & (e_s > 0)))
            d_c = int(np.sum((t_s == ut) & (e_s == target_cause)))
            c_pts.append(c_pts[-1] + S_prev * (d_c / at_risk))
            S_prev *= (1.0 - d_all / at_risk)
            t_pts.append(float(ut))
        t_pts_arr, c_pts_arr = np.array(t_pts), np.array(c_pts)
        idx = int(np.clip(np.searchsorted(t_pts_arr, t_ref, side="right") - 1, 0, len(c_pts_arr) - 1))
        return float(c_pts_arr[idx])

    @staticmethod
    def fit_cox_csh(df_ref: pd.DataFrame, df_arm: pd.DataFrame, target_cause: int = 1, l2_penalty: float = 0.0) -> float:
        df = pd.concat([
            df_ref[["time", "event"]].assign(x=0.0),
            df_arm[["time", "event"]].assign(x=1.0),
        ], ignore_index=True).sort_values("time", ascending=False, kind="mergesort").reset_index(drop=True)
        e = (df["event"].to_numpy() == target_cause).astype(float)
        x = df["x"].to_numpy(dtype=float)
        if e[x == 1.0].sum() < 3:
            return 0.001

        def nll(b_arr):
            b = b_arr[0]
            eb = np.exp(np.clip(b * x, -25.0, 25.0))
            rs = np.cumsum(eb)
            m = e == 1.0
            return -np.sum(b * x[m] - np.log(rs[m])) + l2_penalty * (b**2)

        res = minimize(nll, x0=[0.0], method="BFGS")
        return float(np.exp(res.x[0]))

    @staticmethod
    def _km_censoring_left_limit(times: np.ndarray, events: np.ndarray):
        """Kaplan-Meier estimate of censoring survival G(t-) matching cmprsk::crr."""
        order = np.argsort(times)
        t_o, ev_o = times[order], events[order]
        uniq_c = np.unique(t_o[ev_o == 0])
        g_curr = 1.0
        g_times, g_vals = [0.0], [1.0]
        for ut in uniq_c:
            n_r = np.sum(t_o >= ut)
            d_c = np.sum((t_o == ut) & (ev_o == 0))
            if n_r > 0:
                g_curr *= (1.0 - d_c / n_r)
                g_times.append(float(ut))
                g_vals.append(max(float(g_curr), 1e-8))
        g_times_arr = np.array(g_times, dtype=float)
        g_vals_arr = np.array(g_vals, dtype=float)

        def eval_left(val_arr):
            # Left-continuous limit G(t-) uses side='left' - 1
            idx = np.clip(np.searchsorted(g_times_arr, val_arr, side="left") - 1, 0, len(g_vals_arr) - 1)
            return g_vals_arr[idx]

        return eval_left

    @classmethod
    def fit_fine_gray_shr(
        cls,
        df_ref: pd.DataFrame,
        df_arm: pd.DataFrame,
        target_cause: int = 1,
        l2_penalty: float = 0.0,
    ) -> float:
        """Exact vectorized Fine-Gray subdistribution hazard ratio (IPCW) matching cmprsk::crr."""
        df = pd.concat([
            df_ref[["time", "event"]].assign(x=0.0),
            df_arm[["time", "event"]].assign(x=1.0),
        ], ignore_index=True)
        t = df["time"].to_numpy(dtype=float)
        ev = df["event"].to_numpy(dtype=int)
        x = df["x"].to_numpy(dtype=float)

        if np.sum((ev == target_cause) & (x == 1.0)) < 3:
            return 0.001

        G_left = cls._km_censoring_left_limit(t, ev)

        # Sort by time ascending for O(N log N) exact risk-set sums
        order = np.argsort(t, kind="mergesort")
        t_s, ev_s, x_s = t[order], ev[order], x[order]
        G_s = G_left(t_s)

        ev_mask = ev_s == target_cause
        comp_mask = (ev_s != 0) & (ev_s != target_cause)
        t_ev = t_s[ev_mask]
        x_ev = x_s[ev_mask]
        G_ev = G_s[ev_mask]

        # Precompute exact group counts (x=0 vs x=1) for alive risk sets (t >= t_j)
        # and IPCW-weighted competing risk sets (comp & t < t_j)
        idx_first_ge = np.searchsorted(t_s, t_ev, side="left")
        # Suffix sums of indicators for x==0 and x==1
        is_x0 = (x_s == 0.0).astype(float)
        is_x1 = (x_s == 1.0).astype(float)
        suf_x0 = np.cumsum(is_x0[::-1])[::-1]
        suf_x1 = np.cumsum(is_x1[::-1])[::-1]
        n0_alive = suf_x0[idx_first_ge]
        n1_alive = suf_x1[idx_first_ge]

        # Prefix sums of (1 / G(t_i-)) for competing events strictly before t_j
        invG_comp0 = np.where(comp_mask & (x_s == 0.0), 1.0 / G_s, 0.0)
        invG_comp1 = np.where(comp_mask & (x_s == 1.0), 1.0 / G_s, 0.0)
        pref_c0 = np.concatenate([[0.0], np.cumsum(invG_comp0)])
        pref_c1 = np.concatenate([[0.0], np.cumsum(invG_comp1)])
        w0_comp = G_ev * pref_c0[idx_first_ge]
        w1_comp = G_ev * pref_c1[idx_first_ge]

        R0 = n0_alive + w0_comp
        R1 = n1_alive + w1_comp

        def nll_fg(b_arr):
            b = b_arr[0]
            eb = np.exp(np.clip(b, -25.0, 25.0))
            denom = R0 + R1 * eb
            return -np.sum(b * x_ev - np.log(np.maximum(denom, 1e-12))) + l2_penalty * (b**2)

        def grad_fg(b_arr):
            b = b_arr[0]
            eb = np.exp(np.clip(b, -25.0, 25.0))
            denom = np.maximum(R0 + R1 * eb, 1e-12)
            g = -np.sum(x_ev - (R1 * eb) / denom) + 2.0 * l2_penalty * b
            return np.array([g])

        res = minimize(nll_fg, x0=[0.0], jac=grad_fg, method="BFGS")
        return float(np.exp(res.x[0]))

    def run_causal_audit(self, E_vals: np.ndarray, n_per_arm: int = 2000, t_ref: float = 100.0) -> pd.DataFrame:
        df_h = self.simulate_virtual_cohort(E_vals, n_per_arm=n_per_arm, mode="hormese", seed=42)
        df_m = self.simulate_virtual_cohort(E_vals, n_per_arm=n_per_arm, mode="mascaramento", seed=42)
        df_h0 = df_h[df_h["E"] == 0.0]
        df_m0 = df_m[df_m["E"] == 0.0]
        cif_gen_0 = self.aalen_johansen_cif(df_h0["time"].to_numpy(), df_h0["event"].to_numpy(), 1, t_ref)

        rows = []
        for e_val in E_vals:
            sub_h, sub_m = df_h[df_h["E"] == e_val], df_m[df_m["E"] == e_val]
            cif_g_h = self.aalen_johansen_cif(sub_h["time"].to_numpy(), sub_h["event"].to_numpy(), 1, t_ref)
            cif_a_h = self.aalen_johansen_cif(sub_h["time"].to_numpy(), sub_h["event"].to_numpy(), 2, t_ref)
            cif_g_m = self.aalen_johansen_cif(sub_m["time"].to_numpy(), sub_m["event"].to_numpy(), 1, t_ref)

            hr_csh_h = 1.0 if e_val == 0.0 else self.fit_cox_csh(df_h0, sub_h, 1)
            hr_csh_m = 1.0 if e_val == 0.0 else self.fit_cox_csh(df_m0, sub_m, 1)
            shr_fg_h = 1.0 if e_val == 0.0 else self.fit_fine_gray_shr(df_h0, sub_h, 1)
            shr_fg_m = 1.0 if e_val == 0.0 else self.fit_fine_gray_shr(df_m0, sub_m, 1)

            efs_h = float(np.mean(sub_h["time"] >= t_ref))
            efs_m = float(np.mean(sub_m["time"] >= t_ref))

            red_total = max(cif_gen_0 - cif_g_h, 0.0)
            red_bio = max(cif_g_m - cif_g_h, 0.0)
            phi_hormese = float(np.clip(red_bio / red_total, 0.0, 1.0)) if red_total > 0.02 else 0.0

            rows.append({
                "E": float(e_val),
                "HR_CSH_Horm": hr_csh_h,
                "sHR_FG_Horm": shr_fg_h,
                "HR_CSH_Mask": hr_csh_m,
                "sHR_FG_Mask": shr_fg_m,
                "CIF_gen_Horm": cif_g_h,
                "CIF_amb_Horm": cif_a_h,
                "CIF_gen_Mask": cif_g_m,
                "EFS_Horm": efs_h,
                "EFS_Mask": efs_m,
                "Phi_Antag_Real": phi_hormese,
            })
        return pd.DataFrame(rows)
