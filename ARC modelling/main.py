import pathview
import json
from pathsim.blocks import Scope
import matplotlib.pyplot as plt

# read graph data from a JSON file
with open("arc_same_as_meschini.json", "r") as f:
    graph_data = json.load(f)

model, duration = pathview.make_pathsim_model(graph_data)

model.run(20 * 24 * 3600)

scopes = [block for block in model.blocks if isinstance(block, Scope)]

for i, scope in enumerate(scopes):
    sim_time, data = scope.read()
    plt.figure(i)
    for p, d in enumerate(data):
        lb = scope.labels[p] if p < len(scope.labels) else f"port {p}"
        plt.plot(sim_time / 3600 / 24, d, label=lb)

    plt.xlabel("Time (days)")
    plt.ylabel("inventory (kg)")
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim(bottom=1e-3, top=1)
    plt.xlim(left=1e-1)
plt.show()
