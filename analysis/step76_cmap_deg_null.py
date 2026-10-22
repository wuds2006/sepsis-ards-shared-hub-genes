#!/usr/bin/env python3
"""
Step 76 -- CMap 阴性对照背景 C: DEG-matched null (M4, 严格审稿第 4 轮)

背景 A (非 hub 候选基因) 与真实签名共享 WGCNA 模块结构; 背景 B (全转录组) 则
既无疾病相关性也无签名结构。两者各只提供单侧边界。背景 C 从两病各自上调
DEG (logFC>0, BH FDR<0.05) 的交集中随机抽取 6 基因集 ×100, 同时匹配疾病
相关性与上调签名方向, 用于分离 "六基因特异" 与 "任何协调炎症签名皆然"。

仅查询 Old_CMAP_up (与真实签名的上调查询方向一致)。
输出: 09_docking_drug/cmap_deg_matched_null.csv + summary 行
"""
import csv, json, random, time, urllib.parse, urllib.request

BASE = "/Users/mac/WorkBuddy/2026-08-23-21-32-45/ARDS_Sepsis_reanalysis"
ENRICHR = "https://maayanlab.cloud/Enrichr"
LIB = "Old_CMAP_up"
HUB6 = ["CDK5RAP2", "CD177", "GRB10", "UPP1", "EXOSC4", "VNN1"]
GLUCO = ["dexamethasone", "prednisolone", "prednisone", "hydrocortisone",
         "cortisone", "betamethasone", "triamcinolone", "methylprednisolone",
         "budesonide", "fludrocortisone", "corticosterone"]
REAL_P = 3.66e-4
random.seed(20260917)


def add_list(genes):
    boundary = "----wb76"
    body = ("\r\n".join([
        f"--{boundary}",
        'Content-Disposition: form-data; name="list"',
        "", "\n".join(genes),
        f"--{boundary}",
        'Content-Disposition: form-data; name="description"',
        "", "cmap_deg_matched_null",
        f"--{boundary}--", ""])).encode()
    req = urllib.request.Request(f"{ENRICHR}/addList", data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)["userListId"]


def enrich(uid, lib):
    url = f"{ENRICHR}/enrich?userListId={uid}&backgroundType={urllib.parse.quote(lib)}"
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r).get(lib, [])


def evaluate(terms):
    top10 = terms[:10]
    g_top10 = any(any(g in t[1].lower() for g in GLUCO) for t in top10)
    gluco_rows = [t for t in terms if any(g in t[1].lower() for g in GLUCO)]
    best_gluco_p = min((t[2] for t in gluco_rows), default=None)
    best_p = min((t[2] for t in terms), default=None)
    return g_top10, best_gluco_p, best_p


def load_shared_up():
    def updeg(f):
        ups = set()
        with open(f) as fh:
            for r in csv.DictReader(fh):
                if float(r["logFC"]) > 0 and float(r["adj.P.Val"]) < 0.05:
                    ups.add(r["gene"])
        return ups
    a = updeg(f"{BASE}/01_DEG/GSE76293_DEG.csv")
    b = updeg(f"{BASE}/01_DEG/GSE69528_DEG.csv")
    return sorted((a & b) - set(HUB6))


pool = load_shared_up()
print(f"背景 C 池: {len(pool)} 个共享上调 DEG (去 hub6)", flush=True)

out = f"{BASE}/09_docking_drug/cmap_deg_matched_null.csv"
rows, n_gluco_top10, n_reach = [], 0, 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["set_id", "genes", "gluco_in_top10", "best_gluco_p",
                "best_gluco_reaches_real", "best_p_any_term"])
    for k in range(100):
        genes = random.sample(pool, 6)
        for attempt in range(3):
            try:
                uid = add_list(genes)
                terms = enrich(uid, LIB)
                break
            except Exception as e:
                print(f"  set {k} attempt {attempt}: {e}; retry", flush=True)
                time.sleep(5)
        else:
            raise RuntimeError(f"set {k} failed")
        g10, bgp, bp = evaluate(terms)
        reach = bgp is not None and bgp <= REAL_P
        w.writerow([f"DEG_NULL_{k+1:03d}", ";".join(genes), g10,
                    f"{bgp:.4g}" if bgp is not None else "",
                    reach, f"{bp:.4g}" if bp is not None else ""])
        fh.flush()
        n_gluco_top10 += g10
        n_reach += reach
        if (k + 1) % 10 == 0:
            print(f"  {k+1}/100 done | gluco_top10={n_gluco_top10} | reach_real={n_reach}", flush=True)
        time.sleep(1)

print("\n===== 背景 C (DEG-matched) 结果 =====")
print(f"100 随机集中: top-10 含糖皮质激素 = {n_gluco_top10}/100")
print(f"             达到真实 P ({REAL_P:.2e}) = {n_reach}/100")
with open(f"{BASE}/09_docking_drug/cmap_negative_control_summary.csv", "a", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["background_C_DEG_matched", "Old_CMAP_up", n_reach, n_gluco_top10])
print("summary 已追加")
