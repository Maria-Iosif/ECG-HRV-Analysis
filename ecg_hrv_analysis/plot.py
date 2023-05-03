import matplotlib.pyplot as plt
import seaborn as sns


def r_total(paper, method):
    fig = plt.figure(figsize=(10, 7))

    plt.scatter(x, paper, color="purple", label="total number from paper")
    plt.scatter(x, method, color="purple", alpha=0.5, label="total number from method")
    plt.grid("on")
    plt.xticks(rotation=45)
    plt.legend()
    plt.title("Total number of detected R-peaks")
    plt.ylabel("number of R-peaks")
    plt.xlabel("Record number of subject")


def tp_fp_fn(fpfn_dict):
    avg_dict = avg_fpfn(fpfn_dict)

    data = [avg_dict["tp"], avg_dict["fp"], avg_dict["fn"]]

    fig = plt.figure(figsize=(10, 7))
    x_ax = ["TP", "FP", "FN"]
    ax = sns.boxplot(data, palette="husl")
    ticks = [0, 1, 2]
    ax.set_xticklabels(x_ax)
