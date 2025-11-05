import numpy as np
import matplotlib.pyplot as plt
import bubble_col_glc

def run_scan(diameters, lengths, base_params, elements=50):
    results = {}
    for L in lengths:
        effs = []
        print(f"Scanning L = {L} m")
        for D in diameters:
            params = base_params.copy()
            params.update({"L": L, "D": D, "elements": elements})
            try:
                res = bubble_col_glc.solve(params)
            except Exception as e:
                print(f"  solve failed for D={D:.3f} L={L:.2f}: {e}")
                res = None
            eff = res.get("extraction_efficiency [%]", np.nan) / 100 if res is not None else np.nan
            effs.append(eff)
        results[L] = np.array(effs)
    return results

def plot_results(diameters, results, out_file="extraction_vs_diameter.png"):
    plt.figure(figsize=(8,5))
    for L, effs in results.items():
        plt.plot(diameters, effs, label=f"L = {L} m", lw=1.5)
    plt.xlabel("Column diameter (m)")
    plt.ylabel("Extraction efficiency (fraction)")
    plt.title("Parametric scan: extraction efficiency vs diameter")
    plt.ylim(0, 1)
    plt.xlim(diameters.min(), diameters.max())
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=200)
    print(f"Saved plot to {out_file}")
    try:
        plt.show()
    except Exception:
        pass

def main():
    # base parameters (modify as required)
    base_params = {
        "c_T_inlet": 1.96e-2,   # mol/m^3
        "y_T2_in": 0.0,
        "P_0": 5e5,
        "BCs": "C-C",
        "Flow_l": 560,         # kg/s
        "Flow_g": 10,          # mol/s
        "T": 623,
    }

    diameters = np.linspace(0.05, 4, 40)   # m
    lengths = [1.0, 2.0, 3.0, 4.0]           # m
    elements = 50

    results = run_scan(diameters, lengths, base_params, elements=elements)
    plot_results(diameters, results)

if __name__ == "__main__":
    main()