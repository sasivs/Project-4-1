import os
import pandas as pd
import matplotlib.pyplot as plt
import math as m

def extract_points_for_graph():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\Gplus\Fisher impl\\updated")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\IMDB\Fisher impl\\updated")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\Gplus\Fisher clip\\")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\Gplus")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\Gplus\ep-1\p2-10\\6")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\Gplus\ep-1\p2-NodeNum\\3")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\IMDB\ep-1\p2-NodeNum\\3")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\IMDB\ep-1")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\Gplus\eps_1")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\Gplus\ep-1\p2-10\\6")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\Gplus\ep-1\LocalRRComp\p2-10\\0")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB\ep-1\LocalRRComp\p2-10\\0")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\Gplus\ep-1\p2-NodeNum")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB\ep-1\p2-10\\0_5")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB\ep-1\p2-NodeNum")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\IMDB\eps_1")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\Gplus\eps-1")
    output_direc = os.path.join(direc, "Points")
    # data_noisy = {0:[0]*10, 0.5:[0]*10, 1:[0]*10, 1.5:[0]*10}
    # data_emp = {0:[0]*10, 0.5:[0]*10, 1:[0]*10, 1.5:[0]*10}
    # data_noisy = {0:[0]*4, 0.5:[0]*4, 1:[0]*4, 1.5:[0]*4}
    # data_emp = {0:[0]*4, 0.5:[0]*4, 1:[0]*4, 1.5:[0]*4}
    # data_noisy = {0:[0]*10, 0.5:[0]*10, 1:[0]*10, 1.5:[0]*10, 2.5:[0]*10, 2:[0]*10, 3:[0]*10}
    # data_emp = {0:[0]*10, 0.5:[0]*10, 1:[0]*10, 1.5:[0]*10, 2.5:[0]*10, 2:[0]*10, 3:[0]*10}
    # data_noisy = {3:[0]*10}
    data_samp_emp = {0:[0]*10}
    data_emp1 = {0:[0]*10}
    data_emp2 = {0:[0]*10}
    data_noisy = {0:[0]*10}
    # data_emp = {0:[0]*10}
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi.endswith(".csv"):
                print(fi)
                nodes = int(fi[fi.index("_n")+3:fi.index("_alg")])
                index = int(m.ceil(nodes/1000)) - 1
                # samp_prob = float(fi[fi.index("_p2_")+4:fi.index("_itr")])
                samp_prob = float(fi[fi.index("_p2_")+4:fi.index("_p2flag_")])
                # eps = float(fi[fi.index("_eps_")+5:fi.index("_itr_")])
                # eps = float(fi[fi.index("_eps_")+5:fi.index("_eps_1_0.1")])
                itr = int(fi[fi.index("_itr") + 4:fi.index(".csv")].split("-")[0])-1
                df = pd.read_csv(os.path.join(direc, fi))
                # data_emp[samp_prob][index] = df.loc[itr+2, "Triangle(emp-est)"]
                # data_emp[eps][index] = df.loc[itr+2, "Triangle(emp-est)"]
                # data_emp[eps][index] = df.loc[itr+2, "Trianlges(emp-est)"]
                data_noisy[samp_prob][index] = df.loc[itr+2, "Triangle(est)"]
                # data_noisy[eps][index] = df.loc[itr+2, "Triangle(est)"]
                # data_noisy[eps][index] = df.loc[itr+2, "Trianlges(est)"]
                data_samp_emp[samp_prob][index] = df.loc[itr+2, "Triangle(samp-emp)"]
                data_emp1[samp_prob][index] = df.loc[itr+2, "Triangle(emp1)"]
                data_emp2[samp_prob][index] = df.loc[itr+2, "Triangle(emp2)"]

    df_noisy = pd.DataFrame(data_noisy)
    # df_emp = pd.DataFrame(data_emp)
    df_samp_emp = pd.DataFrame(data_samp_emp)
    df_emp1 = pd.DataFrame(data_emp1)
    df_emp2 = pd.DataFrame(data_emp2)
    df_noisy["Nodes"] = [_*1000 for _ in range(1,11)]
    # df_emp["Nodes"] = [_*1000 for _ in range(1,11)] 
    df_samp_emp["Nodes"] = [_*1000 for _ in range(1,11)]
    df_emp1["Nodes"] = [_*1000 for _ in range(1,11)] 
    df_emp2["Nodes"] = [_*1000 for _ in range(1,11)] 
    # df_noisy["Nodes"] = [_*100000 for _ in range(1,5)]
    # df_emp["Nodes"] = [_*100000 for _ in range(1,5)]
    df_noisy.to_csv(os.path.join(output_direc, "Noisy_points.csv"), index=False)
    # df_emp.to_csv(os.path.join(output_direc, "Emp_points.csv"), index=False)
    df_samp_emp.to_csv(os.path.join(output_direc, "Samp_Emp_points.csv"), index=False)
    df_emp1.to_csv(os.path.join(output_direc, "Emp1_points.csv"), index=False)
    df_emp2.to_csv(os.path.join(output_direc, "Emp2_points.csv"), index=False)
    print(data_noisy)
    # print(data_emp)

# extract_points_for_graph()

def extract_points_for_graph_cliq():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\Gplus\ep-1")
    output_direc = os.path.join(direc, "Points")
    # data_emp3 = {1:[0]*10}
    # data_emp1 = {1:[0]*10}
    # data_emp2 = {1:[0]*10}
    # data_clip3 = {1:[0]*10}
    # data_clip1 = {1:[0]*10}
    # data_clip2 = {1:[0]*10}
    data_emp3 = {1:[0]*11}
    data_emp1 = {1:[0]*11}
    data_emp2 = {1:[0]*11}
    data_clip3 = {1:[0]*11}
    data_clip1 = {1:[0]*11}
    data_clip2 = {1:[0]*11}
    for r,d,f in os.walk(direc):
        for fi in f:
            print(fi)
            if fi.endswith(".csv") and fi.startswith("res"):
                print(fi)
                nodes = int(fi[fi.index("_n")+3:fi.index("_alg")])
                print(nodes)
                index = int(m.ceil(nodes/10000)) - 1
                # eps = float(fi[fi.index("_eps_")+5:fi.index("_itr_")])
                eps = float(fi[fi.index("_eps_")+5:fi.index("_eps_1_0.1")])
                itr = int(fi[fi.index("_itr") + 4:fi.index(".csv")].split("-")[0])-1
                print(itr)
                df = pd.read_csv(os.path.join(direc, fi), index_col=False)
                # df.index = pd.to_numeric(df.index)
                print(df)
                print(df.index.dtype)
                print(df.loc[6])
                data_emp1[eps][index] = df.loc[itr+2, "4-Clique(emp1)"]
                data_emp2[eps][index] = df.loc[itr+2, "4-Clique(emp2)"]
                data_emp3[eps][index] = df.loc[itr+2, "4-Clique(emp3)"]
                data_clip1[eps][index] = df.loc[itr+2, "4-Clique(clip1)"]
                data_clip2[eps][index] = df.loc[itr+2, "4-Clique(clip2)"]
                data_clip3[eps][index] = df.loc[itr+2, "4-Clique(clip3)"]

    df_emp3 = pd.DataFrame(data_emp3)
    df_emp1 = pd.DataFrame(data_emp1)
    df_emp2 = pd.DataFrame(data_emp2)
    df_clip3 = pd.DataFrame(data_clip3)
    df_clip1 = pd.DataFrame(data_clip1)
    df_clip2 = pd.DataFrame(data_clip2)
    print(df_clip3)
    # df_emp1["Nodes"] = [_*10000 for _ in range(1,11)] 
    # df_emp2["Nodes"] = [_*10000 for _ in range(1,11)] 
    # df_emp3["Nodes"] = [_*10000 for _ in range(1,11)] 
    # df_clip1["Nodes"] = [_*10000 for _ in range(1,11)] 
    # df_clip2["Nodes"] = [_*10000 for _ in range(1,11)] 
    # df_clip3["Nodes"] = [_*10000 for _ in range(1,11)] 
    df_emp1["Nodes"] = [_*10000 for _ in range(1,11)] + [107614]
    df_emp2["Nodes"] = [_*10000 for _ in range(1,11)] + [107614]
    df_emp3["Nodes"] = [_*10000 for _ in range(1,11)] + [107614]
    df_clip1["Nodes"] = [_*10000 for _ in range(1,11)] + [107614]
    df_clip2["Nodes"] = [_*10000 for _ in range(1,11)] + [107614]
    df_clip3["Nodes"] = [_*10000 for _ in range(1,11)] + [107614]
    df_emp3.to_csv(os.path.join(output_direc, "Emp3_points.csv"), index=False)
    df_emp1.to_csv(os.path.join(output_direc, "Emp1_points.csv"), index=False)
    df_emp2.to_csv(os.path.join(output_direc, "Emp2_points.csv"), index=False)
    df_clip3.to_csv(os.path.join(output_direc, "clip3_points.csv"), index=False)
    df_clip1.to_csv(os.path.join(output_direc, "clip1_points.csv"), index=False)
    df_clip2.to_csv(os.path.join(output_direc, "clip2_points.csv"), index=False)
    # print(data_noisy)
    # print(data_emp)

# extract_points_for_graph_cliq()

def plot_points():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\IMDB\Fisher impl\\updated\\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\Gplus\Fisher clip\\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\Gplus\Fisher impl\\updated\\Points")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\Gplus\ep-1\p2-10\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\Gplus\\Points")
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi.endswith(".csv"):
                print(fi)
                df = pd.read_csv(os.path.join(direc,fi))
                print(df.columns)
                # plt.show()
                plt.grid()
                df.plot(kind="line", x="Nodes", y='0.0', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_0.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", y='0.5', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_0.5.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", y='1.0', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_1.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", y='1.5', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_1.5.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", y='2.0', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_2.0.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", y='2.5', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_2.5.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", y='3.0', ylabel="Relative Error")
                fname = fi.split(".")[0]+"_sampProb_3.0.jpg"
                plt.savefig(os.path.join(direc,fname))
                df.plot(kind="line", x="Nodes", ylabel="Relative Error")
                # plt.show()
                fname = fi.split(".")[0]+"_sampProb_all.jpg"
                plt.savefig(os.path.join(direc,fname))

# plot_points()

def append_emp_error():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\Gplus\Fisher multi\\varying samp prob(allNodes)")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","BasePermTriImpl\IMDB\Fisher impl")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\Gplus\eps_1")
    output_direc = os.path.join(direc, "updated")
    # print(os.walk(direc).__dict__)
    # print(direc)
    # print(os.getcwd())
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi.endswith(".csv"):
                print(fi)
                itr = int(fi[fi.index("_itr") + 4:fi.index(".csv")].split("-")[0])-1
                # itr = int(fi[fi.index("_itr") + 4:])-1
                df = pd.read_csv(os.path.join(direc, fi))
                print(df)
                print(itr)
                df.columns = ["Trianlges(True)", "Trianlges(est)", "Trianlges(emp-est)", "Trianlges(rel-err(noisy))", "Trianlges(rel-error(emp))"]
                df.loc[itr+1, "Trianlges(est)"] = "AVG(rel-error(noisy))"
                df.loc[itr+1, "Trianlges(emp-est)"] = "AVG(rel-error(emp))"
                df.loc[0:itr, "Trianlges(rel-error(emp))"] = abs(df.loc[0:itr, "Trianlges(True)"].astype(float) - df.loc[0:itr, "Trianlges(emp-est)"].astype(float))/df.loc[0:itr, "Trianlges(True)"].astype(float)
                df.loc[0:itr, "Trianlges(rel-err(noisy))"] = abs(df.loc[0:itr, "Trianlges(True)"].astype(float) - df.loc[0:itr, "Trianlges(est)"].astype(float))/df.loc[0:itr, "Trianlges(True)"].astype(float)
                df.loc[itr+2, "Trianlges(emp-est)"] = df.loc[0:itr, "Trianlges(rel-error(emp))"].mean()
                df.loc[itr+2, "Trianlges(est)"] = df.loc[0:itr, "Trianlges(rel-err(noisy))"].mean()
                df.to_csv(os.path.join(output_direc, fi), index=None)
                # break
# append_emp_error()
                
def consolidate_points():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\Gplus\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\IMDB\ep-1\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB\ep-1\LocalRRComp\p2-10\\0\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\Gplus\ep-1\Points")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\IMDB\ep-1")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\IMDB\eps_1")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB\ep-1")
    output_direc = direc
    dfs = []
    names = []
    for r,d,f in os.walk(direc):
        if r.endswith("Points"):
            print(r)
            for _,d,f in os.walk(r):
                for fi in f:
                    if fi.endswith(".csv"):
                        print(fi)
                        df = pd.read_csv(os.path.join(r,fi))
                        name = fi[:fi.index(".csv")]
                        names.append(name)
                        columns = [col+"_"+name for col in df.columns if col!="Nodes"] + ["Nodes"]
                        df.columns = columns
                        dfs.append(df)
    consolidated_df = dfs[0]
    for i in range(1, len(dfs)):
        consolidated_df = consolidated_df.join(dfs[i].set_index("Nodes"), on=["Nodes"], how="inner", sort="Nodes")
    print(consolidated_df)
    consolidated_df.to_csv(os.path.join(output_direc, "consolidated_points.csv"), index=False)
# consolidate_points()
    
def plot_consolidate_points():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\Gplus\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\IMDB\eps_1\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\Gplus\ep-1")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\Gplus\ep-1")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1\Points")
    output_direc = direc
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi == "consolidated_points.csv":
                df = pd.read_csv(os.path.join(direc,fi))
                columns = list(df.columns)
                for col in columns:
                    if col == "Nodes": continue
                    df.plot(kind="line", x="Nodes", y=col, ylabel="Relative Error", marker='o', mec='r', mfc='r', color='black')
                    plt.grid()
                    fname = col+".jpg"
                    plt.savefig(os.path.join(output_direc, fname))
                columns.remove("Nodes")
                markers = ['o','>','<', '^','v', 's', 'd', 'h','1','2','3','4']
                marker_colors = ['r', 'b', 'g', 'c', 'm', 'y']
                # print(columns)
                fig = df.plot(kind="line", x="Nodes", y=columns, ylabel="Relative Error", color='black')
                plt.grid()
                for i, line in enumerate(fig.get_lines()):
                    line.set_marker(markers[i])
                    line.set_markeredgecolor(marker_colors[i])
                    line.set_markerfacecolor(marker_colors[i])
                    # print(i)
                fig.legend(handles=fig.get_lines())
                # plt.show()
                fname = "consolidated"+".jpg"
                plt.savefig(os.path.join(output_direc, fname))
                # df.plot(kind="line", x="Nodes", y=col, ylabel="Relative Error", marker='o', mec='r', mfcc='r')
# plot_consolidate_points()
                
    
def plot_summary_points():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\Gplus\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\IMDB\eps_1\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\Gplus\ep-1")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results","Gplus")
    output_direc = direc
    print(direc)
    for f in os.listdir(direc):
        # print(f, f.is_file())
        if os.path.isfile(os.path.join(direc,f)) and f.endswith(".csv") and f.startswith("consolidated-RR"):
            print(f)
            df = pd.read_csv(os.path.join(direc,f))
            columns = list(df.columns)
            # columns.pop(2)
            columns.remove("Nodes")
            markers = ['o','>','<', '^','v', 's', 'd', 'h','1','2','3','4']
            marker_colors = ['r', 'b', 'g', 'c', 'm', 'y']
            # print(columns)
            fig = df.plot(kind="line", x="Nodes", y=columns, ylabel="Relative Error", color='black')
            plt.grid()
            for i, line in enumerate(fig.get_lines()):
                line.set_marker(markers[i])
                line.set_markeredgecolor(marker_colors[i])
                line.set_markerfacecolor(marker_colors[i])
                # print(i)
            fig.legend(handles=fig.get_lines())
            # plt.show()
            # fname = f.replace("consolidated", "summary")
            fname = f.replace("consolidated", "RR_v_LS")
            fname = fname.replace("csv", "jpg")
            plt.savefig(os.path.join(output_direc, fname))
# plot_summary_points()
            
def plot_comparison_points():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\Gplus\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\IMDB\eps_1\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\Gplus\ep-1")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB\ep-1")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1\Points")
    output_direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalShuffle\IMDB")
    est = ['Samp_Emp_points']
    samp_prob = ['3n']
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi == "consolidated_points.csv":
                df = pd.read_csv(os.path.join(direc,fi))
                for prob in samp_prob:
                    cols = []
                    for es in est:
                        cols.append(prob+'_'+es)
                    fig = df.plot(kind="line", x="Nodes", y=cols, ylabel="Relative Error", marker='o', mec='r', mfc='r', color='black')
                    plt.grid()
                    markers = ['o','>','<', '^','v', 's', 'd', 'h','1','2','3','4']
                    marker_colors = ['r', 'b', 'g', 'c', 'm', 'y']
                    for i, line in enumerate(fig.get_lines()):
                        line.set_marker(markers[i])
                        line.set_markeredgecolor(marker_colors[i])
                        line.set_markerfacecolor(marker_colors[i])
                        # print(i)
                    fig.legend(handles=fig.get_lines())
                    # plt.show()
                    # break
                    fname = '_vs_'.join(cols)+".jpg"
                    print(fname)
                    plt.show()
                    plt.savefig(os.path.join(output_direc, fname))

                # columns = list(df.columns)
                # for col in columns:
                    # if col == "Nodes": continue
                    # df.plot(kind="line", x="Nodes", y=col, ylabel="Relative Error", marker='o', mec='r', mfc='r', color='black')
                    # plt.grid()
                    # fname = col+".jpg"
                    # plt.savefig(os.path.join(output_direc, fname))
                # columns.remove("Nodes")
                # markers = ['o','>','<', '^','v', 's', 'd', 'h','1','2','3','4']
                # marker_colors = ['r', 'b', 'g', 'c', 'm', 'y']
                # # print(columns)
                # fig = df.plot(kind="line", x="Nodes", y=columns, ylabel="Relative Error", color='black')
                # plt.grid()
                # for i, line in enumerate(fig.get_lines()):
                #     line.set_marker(markers[i])
                #     line.set_markeredgecolor(marker_colors[i])
                #     line.set_markerfacecolor(marker_colors[i])
                #     # print(i)
                # fig.legend(handles=fig.get_lines())
                # # plt.show()
                # fname = "consolidated"+".jpg"
                # plt.savefig(os.path.join(output_direc, fname))
                # df.plot(kind="line", x="Nodes", y=col, ylabel="Relative Error", marker='o', mec='r', mfcc='r')
# plot_comparison_points()

def plot_points_req():
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalRR\Gplus\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\IMDB\eps_1\Points")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","LocalARR\Gplus\ep-1")
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","Gplus")
    # direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1\Points")
    output_direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","Gplus\codaspy-call-for-posters")
    est = ['Samp_Emp_points']
    samp_prob = ['3n']
    fi = "consolidated_p2_3n_ep_1.csv"
    df = pd.read_csv(os.path.join(direc,fi))
    # for prob in samp_prob:
    #     cols = []
    #     for es in est:
    #         cols.append(prob+'_'+es)
    # cols = ["LShuffle*", "LocalRR"]
    cols = ["LShuffle*", "WShuffle*"]
    fig = df.plot(kind="line", x="Nodes", y=cols, ylabel="Relative Error", marker='o', mec='r', mfc='r', color='black')
    plt.grid()
    markers = ['o','>','<', '^','v', 's', 'd', 'h','1','2','3','4']
    marker_colors = ['r', 'b', 'g', 'c', 'm', 'y']
    for i, line in enumerate(fig.get_lines()):
        line.set_marker(markers[i])
        line.set_markeredgecolor(marker_colors[i])
        line.set_markerfacecolor(marker_colors[i])
            # print(i)
        # fig.get_legend().remove()
    fig.legend(handles=fig.get_lines())
        # plt.show()
        # break
        # fname = '_vs_'.join(cols)+".jpg"
    fname = "LShuffle-WShuffle-gplus.jpg"
    print(fname)
    # plt.show()
    plt.savefig(os.path.join(output_direc, fname))

# plot_points_req()  
        
def cliq_recompute_rel_errors():
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1")
    output_direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","ClawShuffle\IMDB\ep-1")
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi.startswith("res_"):
                df = pd.read_csv(os.path.join(direc,fi), index_col=False)
                cols = ["4-Clique(true)","4-Clique(emp1)","4-Clique(emp2)","4-Clique(emp3)","4-Clique(clip1)","4-Clique(clip2)","4-Clique(clip3)"]
                df = df[cols]
                nodes = int(fi.split("_")[2])
                itr = int(fi.split("_")[-1][3:4])-1
                print(fi)
                # return
                df.loc[0:itr, "4-Clique(rel-err(emp1))"] = abs(abs(df.loc[0:itr, "4-Clique(true)"].astype(float))-abs(df.loc[0:itr, "4-Clique(emp1)"].astype(float)))/df.loc[0:itr, "4-Clique(true)"].astype(float)
                df.loc[0:itr, "4-Clique(rel-err(emp2))"] = abs(abs(df.loc[0:itr, "4-Clique(true)"].astype(float))-abs(df.loc[0:itr, "4-Clique(emp2)"].astype(float)))/df.loc[0:itr, "4-Clique(true)"].astype(float)
                df.loc[0:itr, "4-Clique(rel-err(emp3))"] = abs(abs(df.loc[0:itr, "4-Clique(true)"].astype(float))-abs(df.loc[0:itr, "4-Clique(emp3)"].astype(float)))/df.loc[0:itr, "4-Clique(true)"].astype(float)
                df.loc[0:itr, "4-Clique(rel-err(clip1))"] = abs(abs(df.loc[0:itr, "4-Clique(true)"].astype(float))-abs(df.loc[0:itr, "4-Clique(clip1)"].astype(float)))/df.loc[0:itr, "4-Clique(true)"].astype(float)
                df.loc[0:itr, "4-Clique(rel-err(clip2))"] = abs(abs(df.loc[0:itr, "4-Clique(true)"].astype(float))-abs(df.loc[0:itr, "4-Clique(clip2)"].astype(float)))/df.loc[0:itr, "4-Clique(true)"].astype(float)
                df.loc[0:itr, "4-Clique(rel-err(clip3))"] = abs(abs(df.loc[0:itr, "4-Clique(true)"].astype(float))-abs(df.loc[0:itr, "4-Clique(clip3)"].astype(float)))/df.loc[0:itr, "4-Clique(true)"].astype(float)
                
                df.loc[itr+2, "4-Clique(emp1)"] = df.loc[0:itr, "4-Clique(rel-err(emp1))"].mean()
                df.loc[itr+2, "4-Clique(emp2)"] = df.loc[0:itr, "4-Clique(rel-err(emp2))"].mean()
                df.loc[itr+2, "4-Clique(emp3)"] = df.loc[0:itr, "4-Clique(rel-err(emp3))"].mean()
                df.loc[itr+2, "4-Clique(clip1)"] = df.loc[0:itr, "4-Clique(rel-err(clip1))"].mean()
                df.loc[itr+2, "4-Clique(clip2)"] = df.loc[0:itr, "4-Clique(rel-err(clip2))"].mean()
                df.loc[itr+2, "4-Clique(clip3)"] = df.loc[0:itr, "4-Clique(rel-err(clip3))"].mean()
                # print(df.iloc[0:itr+1]["4-Clique(rel-err(emp1))"])
                # print(df.index)
                df.to_csv(os.path.join(output_direc, fi), index=None)
                # break
                # return

def tri_recompute_rel_errors():
    direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\Gplus\eps_1")
    output_direc = os.path.join("D:\Documents\mini project-sem6\ldp implementation\Results\\","WShuffle_star\Gplus\eps_1")
    for r,d,f in os.walk(direc):
        for fi in f:
            if fi.startswith("res_"):
                df = pd.read_csv(os.path.join(direc,fi), index_col=False)
                cols = ["Trianlges(True)","Trianlges(est)","Trianlges(emp-est)","Trianlges(rel-err(noisy))","Trianlges(rel-error(emp))"]
                df = df[cols]
                nodes = int(fi.split("_")[2])
                itr = int(fi.split("_")[-1][3:4])-1
                print(fi)
                # return
                df.loc[0:itr, "Trianlges(rel-err(noisy))"] = abs(abs(df.loc[0:itr, "Trianlges(True)"].astype(float))-abs(df.loc[0:itr, "Trianlges(est)"].astype(float)))/df.loc[0:itr, "Trianlges(True)"].astype(float)
                df.loc[0:itr, "Trianlges(rel-err(emp))"] = abs(abs(df.loc[0:itr, "Trianlges(True)"].astype(float))-abs(df.loc[0:itr, "Trianlges(emp-est)"].astype(float)))/df.loc[0:itr, "Trianlges(True)"].astype(float)
                                
                df.loc[itr+2, "Trianlges(est)"] = df.loc[0:itr, "Trianlges(rel-err(noisy))"].mean()
                df.loc[itr+2, "Trianlges(emp-est)"] = df.loc[0:itr, "Trianlges(rel-err(emp))"].mean()
                # print(df.iloc[0:itr+1]["4-Clique(rel-err(emp1))"])
                # print(df.index)
                df.to_csv(os.path.join(output_direc, fi), index=None)
                # break
                # return

tri_recompute_rel_errors()
    
