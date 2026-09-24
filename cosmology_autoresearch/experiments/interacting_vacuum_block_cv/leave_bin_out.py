#!/usr/bin/env python3
"""Leave-one-redshift-bin-out conditional predictive check for BAO LCDM vs IVS."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular, cho_factor, cho_solve
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "experiments/interacting_vacuum_screen"))
import background_bao as bao
import interacting_vacuum_profile as iv

TOL = 1e-8


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load():
    return bao.BAOData.from_files(ROOT/"context/data/desi_dr2_mean.txt", ROOT/"context/data/desi_dr2_cov.txt")


def subset(data, ix):
    ix = np.asarray(ix, dtype=int)
    return bao.BAOData(data.z[ix], data.value[ix], data.observable[ix], data.covariance[np.ix_(ix, ix)])


def project_alpha(data, chol, q):
    yw = solve_triangular(chol, data.value, lower=True, check_finite=False)
    qw = solve_triangular(chol, q, lower=True, check_finite=False)
    alpha = float(np.clip(np.dot(qw, yw)/np.dot(qw, qw), *iv.ALPHA_BOUNDS))
    r = yw-alpha*qw
    return float(r@r), alpha, alpha*q


def predict_lcdm(data, om):
    return bao.predict_bao(data.z, data.observable, 1.0, float(om))


def full_domain_prediction(data, om, g):
    # Appending z=2.33 forces integrate_profile to validate the whole required domain,
    # including folds which have omitted the highest-redshift data group.
    zall = np.concatenate((data.z, np.array([0.0, iv.ZMAX])))
    bg = iv.integrate_profile(om, g, zall)
    dm, e = bg["distance"][:len(data.z)], bg["e"][:len(data.z)]
    dh = 1/e
    dv = np.cbrt(data.z*dm*dm*dh)
    q = np.empty_like(data.value, dtype=np.float64)
    for label, vals in (("DM_over_rs", dm), ("DH_over_rs", dh), ("DV_over_rs", dv)):
        q[data.observable == label] = vals[data.observable == label]
    return q, bg


def fit_lcdm(train, starts=13):
    chol = np.linalg.cholesky(train.covariance)
    attempts = []
    for i, om0 in enumerate(np.linspace(*bao.FIT_BOUNDS["lcdm"][1], starts)):
        def obj(x):
            try:
                return project_alpha(train, chol, predict_lcdm(train, x[0]))[0]
            except (ValueError, FloatingPointError, OverflowError):
                return 1e80
        res = minimize(obj, [om0], method="L-BFGS-B", bounds=[bao.FIT_BOUNDS["lcdm"][1]],
                       options={"maxiter":700,"ftol":1e-14,"gtol":1e-9,"maxls":40})
        score, alpha, pred = project_alpha(train, chol, predict_lcdm(train, res.x[0]))
        om = float(res.x[0])
        attempts.append({"start_index":i,"x0":[float(om0)],"x":[om],"chi2_train":score,
                         "alpha":alpha,"success":bool(res.success),"status":int(res.status),
                         "message":str(res.message),"nit":int(res.nit),"nfev":int(res.nfev),
                         "parameter_bound_hits":(["Omega_m:lower"] if abs(om-.05)<TOL else [])+
                             (["Omega_m:upper"] if abs(om-.6)<TOL else [])+
                             (["alpha:lower"] if abs(alpha-iv.ALPHA_BOUNDS[0])<1e-12 else [])+
                             (["alpha:upper"] if abs(alpha-iv.ALPHA_BOUNDS[1])<1e-9 else [])})
    successful=[a for a in attempts if a["success"]]
    best=min(successful or attempts,key=lambda a:a["chi2_train"])
    return {"Omega_m":best["x"][0],"Gamma_over_H0":0.0,"alpha":best["alpha"],
            "chi2_train":best["chi2_train"],"prediction_train":(best["alpha"]*predict_lcdm(train,best["x"][0])).tolist(),
            "attempts":attempts,"selected_success":best["success"],"selected_status":best["status"],
            "selected_start_index":best["start_index"],"parameter_bound_hits":best["parameter_bound_hits"]}


def fit_ivs(train, multistarts=8, grid_n=7):
    chol=np.linalg.cholesky(train.covariance)
    yw=solve_triangular(chol,train.value,lower=True,check_finite=False)
    failures={}
    def score(om,g, keep=False):
        try:
            q,bg=full_domain_prediction(train,float(om),float(g))
            qw=solve_triangular(chol,q,lower=True,check_finite=False)
            alpha=float(np.clip(qw@yw/(qw@qw),*iv.ALPHA_BOUNDS))
            rr=yw-alpha*qw
            ans=float(rr@rr)
            if not np.isfinite(ans): raise ValueError("nonfinite profiled score")
            return (ans,alpha,q,bg) if keep else ans
        except (iv.DomainError,ValueError,FloatingPointError,OverflowError) as ex:
            failures[str(ex)]=failures.get(str(ex),0)+1
            return (1e80,None,None,None) if keep else 1e80
    bounds=[iv.OM_BOUNDS,iv.GAMMA_BOUNDS]
    starts=[np.array([om,g]) for om in np.linspace(*iv.OM_BOUNDS,4) for g in np.linspace(*iv.GAMMA_BOUNDS,3)]
    starts += [np.array([.12,-.35]),np.array([.22,.35]),np.array([.42,-.35]),np.array([.52,.35])]
    starts=starts[:multistarts]
    attempts=[]
    for i,x0 in enumerate(starts):
        res=minimize(lambda x:score(x[0],x[1]),x0,method="L-BFGS-B",bounds=bounds,
                     options={"maxiter":1000,"ftol":1e-14,"gtol":2e-8,"maxls":40})
        val,alpha,q,bg=score(*res.x,keep=True)
        om,g=map(float,res.x)
        hits=((["Omega_m:lower"] if abs(om-iv.OM_BOUNDS[0])<TOL else [])+
              (["Omega_m:upper"] if abs(om-iv.OM_BOUNDS[1])<TOL else [])+
              (["Gamma_over_H0:lower"] if abs(g-iv.GAMMA_BOUNDS[0])<TOL else [])+
              (["Gamma_over_H0:upper"] if abs(g-iv.GAMMA_BOUNDS[1])<TOL else [])+
              (["alpha:lower"] if alpha is not None and abs(alpha-iv.ALPHA_BOUNDS[0])<1e-12 else [])+
              (["alpha:upper"] if alpha is not None and abs(alpha-iv.ALPHA_BOUNDS[1])<1e-9 else []))
        attempts.append({"start_index":i,"x0":x0.tolist(),"x":[om,g],"chi2_train":val,
                         "alpha":alpha,"endpoint_valid_physical":alpha is not None,
                         "domain_minima":None if bg is None else {k:bg[k] for k in ("minimum_e2","minimum_matter","minimum_vacuum")},
                         "success":bool(res.success),"status":int(res.status),"message":str(res.message),
                         "nit":int(res.nit),"nfev":int(res.nfev),"parameter_bound_hits":hits})
    # Independent coarse grid/local refinement check for the two shape coordinates.
    omgrid=np.linspace(*iv.OM_BOUNDS,grid_n); ggrid=np.linspace(*iv.GAMMA_BOUNDS,grid_n)
    grid=[]
    for om in omgrid:
        for g in ggrid:
            grid.append((score(om,g),float(om),float(g)))
    top=sorted((x for x in grid if x[0]<1e80),key=lambda x:x[0])[:3]
    refinements=[]
    for j,(_,om0,g0) in enumerate(top):
        res=minimize(lambda x:score(x[0],x[1]),[om0,g0],method="Powell",bounds=bounds,
                     options={"maxiter":700,"xtol":1e-9,"ftol":1e-12})
        val,alpha,q,bg=score(*res.x,keep=True)
        refinements.append({"start_index":j,"grid_start":[om0,g0],"x":[float(x) for x in res.x],
                            "chi2_train":val,"alpha":alpha,"success":bool(res.success),
                            "status":int(res.status),"message":str(res.message),"nit":int(res.nit),"nfev":int(res.nfev),
                            "parameter_bound_hits":(((["Omega_m:lower"] if abs(res.x[0]-iv.OM_BOUNDS[0])<TOL else [])+
                              (["Omega_m:upper"] if abs(res.x[0]-iv.OM_BOUNDS[1])<TOL else [])+
                              (["Gamma_over_H0:lower"] if abs(res.x[1]-iv.GAMMA_BOUNDS[0])<TOL else [])+
                              (["Gamma_over_H0:upper"] if abs(res.x[1]-iv.GAMMA_BOUNDS[1])<TOL else [])))})
    candidates=[a for a in attempts if a["endpoint_valid_physical"]]+[
        {**a,"x":a["x"],"chi2_train":a["chi2_train"],"alpha":a["alpha"],"success":a["success"],
         "status":a["status"],"parameter_bound_hits":a["parameter_bound_hits"]} for a in refinements if a["alpha"] is not None]
    successes=[a for a in candidates if a["success"]]
    best=min(successes or candidates,key=lambda a:a["chi2_train"])
    om,g=best["x"]
    val,alpha,q,bg=score(om,g,keep=True)
    return {"Omega_m":float(om),"Gamma_over_H0":float(g),"alpha":float(alpha),"chi2_train":float(val),
            "prediction_train":(alpha*q).tolist(),"selected_success":best["success"],"selected_status":best["status"],
            "parameter_bound_hits":best["parameter_bound_hits"],"multistart_attempts":attempts,
            "grid_refinement_attempts":refinements,"coarse_grid":{"shape":[grid_n,grid_n],
            "finite_physical_point_count":sum(x[0]<1e80 for x in grid),"best_points":[{"chi2":x[0],"Omega_m":x[1],"Gamma_over_H0":x[2]} for x in top]},
            "domain_failure_evaluations_by_reason":failures,
            "domain_minima":{"minimum_e2":bg["minimum_e2"],"minimum_matter":bg["minimum_matter"],
                             "minimum_vacuum":bg["minimum_vacuum"]}}


def conditional_predict(data, train_idx, held_idx, train_prediction, held_prediction):
    c=data.covariance
    ctt=c[np.ix_(train_idx,train_idx)]
    cht=c[np.ix_(held_idx,train_idx)]
    chh=c[np.ix_(held_idx,held_idx)]
    fac=cho_factor(ctt,lower=True,check_finite=False)
    adjustment=cht@cho_solve(fac,data.value[train_idx]-train_prediction,check_finite=False)
    mean=held_prediction+adjustment
    cov=chh-cht@cho_solve(fac,cht.T,check_finite=False)
    cov=(cov+cov.T)/2
    pred_resid=data.value[held_idx]-mean
    chol=np.linalg.cholesky(cov)
    white=solve_triangular(chol,pred_resid,lower=True,check_finite=False)
    chi2=float(white@white)
    # Independent dense-solve check of both Schur-complement expressions.
    solved_resid=np.linalg.solve(ctt,data.value[train_idx]-train_prediction)
    solved_cross=np.linalg.solve(ctt,cht.T)
    mean_direct=held_prediction+cht@solved_resid
    cov_direct=chh-cht@solved_cross
    direct_chi2=float(pred_resid@np.linalg.solve(cov,pred_resid))
    checks={"max_abs_conditional_mean_difference_vs_dense_solve":float(np.max(np.abs(mean-mean_direct))),
            "max_abs_schur_covariance_difference_vs_dense_solve":float(np.max(np.abs(cov-cov_direct))),
            "conditional_covariance_symmetry_max_abs":float(np.max(np.abs(cov-cov.T))),
            "conditional_covariance_min_eigenvalue":float(np.min(np.linalg.eigvalsh(cov))),
            "chi2_difference_vs_dense_solve":float(abs(chi2-direct_chi2))}
    return mean,cov,chi2,float(np.linalg.slogdet(cov)[1]),checks


def main():
    t0=time.monotonic(); data=load(); folds=[]
    for z in np.unique(data.z):
        held=np.flatnonzero(data.z==z); train=np.flatnonzero(data.z!=z)
        tr=subset(data,train); h=subset(data,held)
        lm=fit_lcdm(tr)
        ivs=fit_ivs(tr)
        lpred=np.asarray(lm["prediction_train"]); ipred=np.asarray(ivs["prediction_train"])
        lhm=bao.predict_bao(h.z,h.observable,lm["alpha"],lm["Omega_m"])
        ihq,_=full_domain_prediction(h,ivs["Omega_m"],ivs["Gamma_over_H0"])
        ihm=ivs["alpha"]*ihq
        lmean,lcov,lchi2,llogdet,lchecks=conditional_predict(data,train,held,lpred,lhm)
        imean,icov,ichi2,ilogdet,ichecks=conditional_predict(data,train,held,ipred,ihm)
        folds.append({"held_out_z":float(z),"held_out_row_indices_0based":held.tolist(),
            "held_out_rows":[{"row_0based":int(i),"z":float(data.z[i]),"observable":str(data.observable[i]),"value":float(data.value[i])} for i in held],
            "training_row_indices_0based":train.tolist(),"lcdm_fit":lm,"ivs_fit":ivs,
            "lcdm_conditional_mean":lmean.tolist(),"ivs_conditional_mean":imean.tolist(),
            "lcdm_conditional_covariance":lcov.tolist(),"ivs_conditional_covariance":icov.tolist(),
            "lcdm_conditional_chi2":lchi2,"ivs_conditional_chi2":ichi2,
            "delta_ivs_minus_lcdm_conditional_chi2":ichi2-lchi2,
            "lcdm_conditional_logdet":llogdet,"ivs_conditional_logdet":ilogdet,
            "lcdm_conditional_numerical_checks":lchecks,"ivs_conditional_numerical_checks":ichecks})
        print(f"z={z:g}: held rows={len(held)} conditional chi2 LCDM={lchi2:.6f} IVS={ichi2:.6f} delta={ichi2-lchi2:+.6f}",flush=True)
    result={"status":"complete_leave_one_redshift_bin_out_conditional_prediction",
      "data_contract":{"row_order":"released desi_dr2_mean.txt line order, 0-based indices","row_count":len(data.z),
        "redshift_groups":[float(x) for x in np.unique(data.z)],"full_released_covariance_used":True,
        "conditional_formula":"mu(H|T)=mu_H+C_HT C_TT^{-1}(y_T-mu_T); S=C_HH-C_HT C_TT^{-1}C_TH; chi2=(y_H-mu(H|T))^T S^{-1}(y_H-mu(H|T))",
        "group_delta_definition":"conditional chi2(IVS)-conditional chi2(flat LCDM); negative favors IVS for that plug-in score",
        "model_fitting":"alpha, Omega_m, and IVS Gamma/H0 fit to training rows only; alpha analytically profiled with training covariance",
        "physical_validation":"every IVS objective candidate checked densely over 0<=z<=2.33 via full_domain_prediction",
        "bounds_are_search_box_not_prior":{"alpha":iv.ALPHA_BOUNDS,"Omega_m":iv.OM_BOUNDS,"Gamma_over_H0":iv.GAMMA_BOUNDS}},
      "folds":folds,"aggregate":{"lcdm_sum_conditional_chi2":float(sum(f["lcdm_conditional_chi2"] for f in folds)),
        "ivs_sum_conditional_chi2":float(sum(f["ivs_conditional_chi2"] for f in folds)),
        "sum_delta_ivs_minus_lcdm":float(sum(f["delta_ivs_minus_lcdm_conditional_chi2"] for f in folds)),
        "interpretation":"sum across seven overlapping leave-group-out predictive assessments; not a joint likelihood or evidence"},
      "limitations":["plug-in parameter fits, not posterior predictive cross-validation","only seven BAO redshift groups",
        "held-out groups share the same survey and released covariance; not an independent survey","late-time compressed BAO screen only; no CMB, sound-horizon calibration, or perturbation test"],
      "runtime_seconds":time.monotonic()-t0,"environment":{"python":sys.version,"numpy":np.__version__,
        "platform":platform.platform(),"threads":{"OMP_NUM_THREADS":os.getenv("OMP_NUM_THREADS"),
        "OPENBLAS_NUM_THREADS":os.getenv("OPENBLAS_NUM_THREADS"),"MKL_NUM_THREADS":os.getenv("MKL_NUM_THREADS"),
        "NUMEXPR_NUM_THREADS":os.getenv("NUMEXPR_NUM_THREADS")}},
      "sha256_inputs":{"context/data/desi_dr2_mean.txt":sha(ROOT/"context/data/desi_dr2_mean.txt"),
        "context/data/desi_dr2_cov.txt":sha(ROOT/"context/data/desi_dr2_cov.txt"),
        "scripts/background_bao.py":sha(ROOT/"scripts/background_bao.py"),
        "experiments/interacting_vacuum_screen/interacting_vacuum_profile.py":sha(ROOT/"experiments/interacting_vacuum_screen/interacting_vacuum_profile.py"),
        "experiments/interacting_vacuum_block_cv/leave_bin_out.py":sha(Path(__file__).resolve())}}
    out=ROOT/"experiments/interacting_vacuum_block_cv/result.json"
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(f"aggregate: LCDM={result['aggregate']['lcdm_sum_conditional_chi2']:.6f} IVS={result['aggregate']['ivs_sum_conditional_chi2']:.6f} delta={result['aggregate']['sum_delta_ivs_minus_lcdm']:+.6f}; runtime={result['runtime_seconds']:.1f}s")


if __name__=="__main__": main()
