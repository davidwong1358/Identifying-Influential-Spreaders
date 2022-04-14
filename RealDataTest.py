import snap
import pandas as pd
import random
import collections

import matplotlib.pyplot as plt

import concurrent.futures
import time

from itertools import combinations

#########################################################
#Project Algorithm

#K-shell Decomposition: return node ks and top n influential spreaders
def kshellDecomp(graph_file, num_initSP):
    graph = snap.LoadEdgeList(snap.TUNGraph, graph_file, 0, 1)
    cp = snap.ConvertGraph(type(graph), graph) #Copy the graph
    degree = 1 #Asumption: All node has been connected, all node have degree >= 1
    
    tmp = []
    shell = []
    
    while True:
        empty = True #Check whether the updated graph contain node with certain degree
        for NI in cp.Nodes():
            #Check if the graph contains node with that degree
            if(NI.GetDeg() <= degree): 
                empty = False
                break

        #No more node with that degree in the graph/ Cannot further delete node
        if empty == True: 
            degree += 1
            shell.append(tmp)
            tmp = []

        if empty == False:
            node_set = []
            for NI in cp.Nodes():
                if(NI.GetDeg() <= degree): #Record the node id which fulfill condition
                    node_set.append(NI.GetId())

            for nd in node_set: #Delete the recorded node and record the id for the shell
                cp.DelNode(nd)
                tmp.append(nd)

        if (cp.GetNodes() == 0): #Final case: All node deleted
            shell.append(tmp)
            break

    #Create dataframe for HGD Algo and obtaining top n influential spreader
    data = []
    ks = 1
    for s in shell:
        for n in s:
            NI = graph.GetNI(n)
            data.append([n, ks, NI.GetDeg()])
        ks += 1

    df = pd.DataFrame(data, columns=['node', 'ks', 'degree']) 
    df = df.set_index('node')
    df.sort_values(by=['node', 'ks', 'degree'] ,inplace=True, ascending=[True, False, False])

    #Obtaining top n influential spreader by greatest ks > degree > node number
    df_cp = df.copy()
    df_cp.sort_values(by=['ks', 'degree', 'node'] ,inplace=True, ascending=[False, False, True])
    
    #Demo: Print result
    print(str(shell))
    max_ks_grp = df_cp[df_cp['ks'] == df_cp['ks'].max()] #Largest ks with sorted by degree
    print(max_ks_grp)

    return df, df_cp.head(num_initSP).index.values.tolist()

#Find the initial spreaders by HGD Algorithm: return top n spreaders
def hgdAlgo(graph_file, df, num_initSP):

    graph = snap.LoadEdgeList(snap.TUNGraph, graph_file, 0, 1)
    df_cp = df.copy()
    grp = []

    #While loop until all node with ks = 0
    while True:
        #Step 1: Find initial group (Initial Spreader + Neighbour with same ks)
        core = []
        tmp = []

        #Find max ks grp by using dataframe
        max_ks_grp = df_cp[df_cp['ks'] == df_cp['ks'].max()] 

        #Return first index with greatest degree in max ks group
        max_dg_node = max_ks_grp[max_ks_grp['degree'] == max_ks_grp['degree'].max()].index.values.tolist()[0]

        NI = graph.GetNI(max_dg_node)
        core.append(max_dg_node)
        df_cp['ks'][max_dg_node] = 0

        #Add neighbour with same ks value
        nlist = [e for e in NI.GetOutEdges()]
        tmp = [item for item in nlist if item in max_ks_grp.index.values.tolist()]
        for n in tmp:
            core.append(n)
            df_cp['ks'][n] = 0

        tmp = []

        #Step 2: While loop until no neighbours can be further added into group
        while True:
            c = core.copy() #Copy for compare the group identity

            #Neighbour of node in initial core
            for n in core:
                #Create list for storing neighbour node
                NI = graph.GetNI(n)
                olist = [e for e in NI.GetOutEdges()] #Out edge = In edge for undirected graph            
                #Add the neighbour nodes for each nodes in core
                for nd in olist:
                    if (df_cp['ks'][nd] != 0 and nd not in tmp):
                        tmp.append(nd)

            #Sort the node by degree in ascending order
            df_filter = df_cp[df_cp.index.isin(tmp)].copy()
            df_filter.sort_values(by=['degree', 'ks', 'node'] , inplace=True, ascending=[True,False,True])
            tmp = df_filter.index.values.tolist()

            #Find connection (degree) between node & group (kin & kout)
            while True:
                nd_set = [] #Store the node that kin >= kout for a iteration
                for i in tmp:
                    NI = graph.GetNI(i)
                    #Neighbour list, common node, different node
                    nlist = [e for e in NI.GetOutEdges()]
                    kin = len([nd for nd in nlist if nd in core]) #Common Node of neighbour node list & core node list (N n C)
                    kout = len(nlist) - kin #Difference of neighbour list & core list (N \ C)

                    if (kin >= kout):
                        nd_set.append(i)
                
                for nd in nd_set:
                    if nd in tmp:
                        tmp.remove(nd)
                        core.append(nd)
                        df_cp['ks'][nd] = 0

                if(len(nd_set) == 0): #Check if no neighbour nodes can be added for a node
                    break

            #Check if no node is further added to core
            if (c == core):
                grp.append(core)
                break
        ############################################

        #Step 4: Check if all node has ks = 0, end the loop if all ks = 0
        if len(df_cp[df_cp['ks'] > 0 ].index) == 0:
            break
       
    #Step 4: Rank the group by node number, and find the initial spreader by top degree of node
    grp.sort(key=len, reverse = True)
    
    init_sp = []

    for c in grp:
        #Filter the row that index is in the core
        df_filter = df_cp[df_cp.index.isin(c)].copy()
        #Get the greatest degree node from filtered row
        sp = df_filter[df_filter['degree'] == df_filter['degree'].max()].index.values.tolist()[0]
        init_sp.append(sp)
    
    #Used for demo
    print(str(grp)) #Test the group

    return init_sp[:num_initSP]

#Simulation Model for testing, i = contact rate, r = recovery rate: return filename, infect node, running time
def sirModel(graph_file, init_SP, time_sir, i, r, output_file):
    print("SIR start: " + str(output_file))
    graph = snap.LoadEdgeList(snap.TUNGraph, graph_file, 0, 1)

    data = [] #Store graph data: node, source node, status (S/I/R)
    infected = [] #Store infected node
    result = [] #Store the result for single simulation

    for NI in graph.Nodes():
        data.append([NI.GetId(), "Susceptible"])

    #Dataframe of node set at latest time
    df = pd.DataFrame(data, columns=['node', 'status']) 
    df = df.set_index('node')

    for n in init_SP: #Add source nodes as infected
        df.loc[n, 'status'] = "Infected"
        infected.append(n)

    #Count the number for S, I, R status for the graph at time = 0
    sus = len(df[df['status'] == "Susceptible"].index)
    inf = len(df[df['status'] == "Infected"].index)
    rec = len(df[df['status'] == "Recovered"].index)
    
    #Time, no of sus, inf, rec nodes
    result.append([0, sus, inf, rec])
    #print("Infection Rate: %f" % (i)) #Debug

    for t in range (1, time_sir + 1):
        initSP_list = infected.copy() #Record the infected node at last time
        random.shuffle(initSP_list) #Shuffle list to avoid head node in infected list always having priority to influence other node
        #print("\nTime: %i, Previous Infected Node: %s" % (t, str(initSP_list))) #Debug

        #Neighbour Infection
        #print("Infection process in time %i" % (t)) #Debug
        for n in initSP_list:
            NI = graph.GetNI(n)
            nlist = [e for e in NI.GetOutEdges() if df['status'][e] == "Susceptible"] #Get the neighbour list that is susceptible
            for e in nlist: 
                p_i = random.uniform(0.00, 1.00) #Generate random number for simulating infection 
                #print("Susceptible node: %s, Infect result: %f" % (str(e), p_i)) #Debug
                if p_i <= i:
                    infected.append(e)
                    #print("Infected: %s" % (str(e)))
                    df.loc[e, 'status'] = "Infected"

        #print("Recovery process in time %i" % (t))
        for n in initSP_list:
            #Recover the node
            if df['status'][n] == "Infected":
                p_r = random.uniform(0.00, 1.00)
                #print("Infected node: %s, Recover result: %f" % (str(n), p_r)) #Debug
                if p_r <= r:
                    infected.remove(n)
                    df.loc[n, 'status'] = "Recovered"
                    #print("Recovered: %s" % (str(n)))

        #Record the no of the status of the node for each time
        sus = len(df[df['status'] == "Susceptible"].index)
        inf = len(df[df['status'] == "Infected"].index)
        rec = len(df[df['status'] == "Recovered"].index)
        result.append([t, sus, inf, rec])
    
    final_inf = len(df[df['status'] == "Infected"].index) + len(df[df['status'] == "Recovered"].index)
    #Display Plotted Graph for SIR
    result_sir = pd.DataFrame(result, columns=['time', 'Susceptible', 'Infected', "Recovered"])
    result_sir.plot(x="time", y=["Susceptible", "Infected", "Recovered"], kind="line")

    #result_sir.to_csv(str(output_file) + '_SIR.csv') #Generate the output
    #plt.savefig(str(output_file) + '_SIR.png') #Store the output graph

    print("SIR end: " + str(output_file))
    return str(output_file), final_inf

#SIR Model for GA (Simplified Version): return infect node, src contribution, tested node set
def simpleSIR(graph_file, init_SP, i, time_sir):
    graph = snap.LoadEdgeList(snap.TUNGraph, graph_file, 0, 1)

    data = [] #Store graph data: node, source node, status (S/I/R)
    infected = []
    src_influ = [] #Store the node influence by src node (Either directly or indirectly)


    for NI in graph.Nodes():
        data.append([NI.GetId(), -1, "Susceptible"])
    
    df = pd.DataFrame(data, columns=['node', 'source', 'status']) 
    df = df.set_index('node')

    #Add source nodes as infected
    for n in init_SP: 
        df.loc[n, 'status'] = "Infected"
        df.loc[n, 'source'] = n
        infected.append(n)

    for t in range (0, time_sir):
        initSP_list = infected.copy() #Record the infected node at last time
        random.shuffle(initSP_list) 

        #Neighbour Infection
        for n in initSP_list:
                NI = graph.GetNI(n)
                nlist = [e for e in NI.GetOutEdges() if df['status'][e] == "Susceptible"] #Get the neighbour list that is susceptible
                for e in nlist: 
                    p_i = random.uniform(0.00, 1.00) #Generate random number for simulating infection 
                    if p_i <= i:
                        infected.append(e)
                        df.loc[e, 'status'] = "Infected"
                        df.loc[e, 'source'] = df.loc[n, 'source']

                #Recover the node
                if df['status'][n] == "Infected":
                        infected.remove(n)
                        df.loc[n, 'status'] = "Recovered"

    for n in init_SP:
        influ = len(df[df['source'] == n].index)
        src_influ.append([n, influ])

    total_inf = len(df[df['status'] == "Infected"].index) + len(df[df['status'] == "Recovered"].index)

    #Return final Recovered node (Infected before) and source node influence
    return total_inf, src_influ, init_SP

#GA (popu = no of soln, chromo: no of node in a chromosome), return fittest node set, fitness, filename, time
def GeneAlgo(graph_file, inf, popu, chromo, max_gen, output_file, time_sir):
    graph = snap.LoadEdgeList(snap.TUNGraph, graph_file, 0, 1)

    print("GA start: " + str(output_file))
    start_time = time.time()
    
    nodes = graph.GetNodes()
    
    generation = 1
    used_set = [] #Record used set to prevent used node set is used in afterward generation

    population = []

    fit_max = -1 #Max fitness in a generation
    fit_set = [] #Max fitness node set in a generation

    chrm_list = [] #Store each chromosome and its fitness for each SIR result
    fit_list = [] #Store fitness for each SIR result
    src_contri = {} #Store source node contribution for a node set (Dict)

    ga_stat = []
    
    #print("Generation: %i" % (generation)) #Debug

    #First Generation: Generate Chromosome by random
    while len(population) < popu:
        chromosome = random.sample(range(0, nodes - 1), chromo)
        chromosome.sort() #Sort the node set for comparing duplicate node
        if chromosome not in used_set:
            population.append(chromosome)
            used_set.append(chromosome)

    #Evaluation of Fitness
    #Run SIR for each chromosome, record fitness of each chromosome and indiv influ for each gene
    #Multithread
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(simpleSIR, graph_file, chrm, inf, time_sir) for chrm in population]
        for future in concurrent.futures.as_completed(futures):
            fitness, src, chrm_1 = future.result()
            #print(str(fitness)+ " " + str(chrm_1))
            chrm_list.append([chrm_1, fitness])
            fit_list.append(fitness)

            #Find best fit node set in one generation
            if fitness > fit_max:
                fit_set = chrm_1
                fit_max = fitness
        
            #Add individual source contribution
            #Use dict, breakdown node as key, src influ as value. Append if no key exist or same key with higher influ
            for n in src:
                key, value = n[0], n[1]
                if key in src_contri:
                    if value >= src_contri[key]:
                        src_contri[key] = value
                else:
                    src_contri[key] = value


    #Rank the source node influrence
    src_contri = dict(sorted(src_contri.items(), key=lambda item: item[1], reverse=True))
    #print("Source Contribution:\n%s" % (str(src_contri))) #Debug
    #Record stat for current generation
    fit_avg = round(sum(fit_list) / len(fit_list), 2)
    ga_stat.append([generation, fit_max, fit_avg, fit_set])

    #Second generation and after
    while generation < max_gen:
        print("End of generation : %d" % (generation))
        generation += 1
        print("\nGeneration: %i" % (generation))
        #Reproduction
        #Create chromosome for next generation
        population = [] #Clear the population list for allocating node set in next generation

        #Best fit set
        population.append(fit_set) #Best fitness set
        used_set.append(fit_set)
        #print("Best Fitness Set: %s" % (str(fit_set))) #Debug

        #Top n influential source
        gene = set()
        for n in range(chromo):
            nd = list(src_contri.keys())[n]
            gene.add(nd)
        gene = list(gene)
        #print("Top Contribution Set: %s" % (str(gene))) #Debug
        gene.sort()
        population.append(gene)
        used_set.append(gene)
        
        

        #Mutation
        #First half: Random in Top half node, Remaining half: Random without duplicate
        top_half = len(src_contri) // 2

        #If population after filled mutated chromosome is odd, filled chromosome pair in crossover step may exceed the population
        mut_half = (popu - 2) // 2
        len_mut = mut_half + 2
        if (popu - len_mut) % 2 == 1:
            len_mut += 1

        while len(population) < popu:
            gene = set()
            while len(gene) < chromo // 2:
                nd = list(src_contri.keys())[random.randrange(top_half)]
                gene.add(nd)
            while len(gene) < chromo:
                gene.add(random.randrange(nodes - 1))

            gene = list(gene) #Debug
            #print ("New Gene : " + str(gene))
            gene.sort()

            if gene not in used_set:
                population.append(gene)
                used_set.append(gene)


        #Reset for current generation
        fit_list = []
        src_contri = {}
        chrm_list = []

        #Evaluation of fitness by Multithread
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(simpleSIR, graph_file, chrm, inf, time_sir) for chrm in population]
            for future in concurrent.futures.as_completed(futures):
                fitness, src, chrm_1 = future.result()
                #print(str(fitness)+ " " + str(chrm_1))
                chrm_list.append([chrm_1, fitness])
                fit_list.append(fitness)

                #Find best fit node set in one generation
                if fitness > fit_max:
                    fit_set = chrm_1
                    fit_max = fitness
        
                #Add individual source contribution
                #Use dict, breakdown node as key, src influ as value. Append if no key exist or same key with higher influ
                for n in src:
                    key, value = n[0], n[1]
                    if key in src_contri:
                        if value >= src_contri[key]:
                            src_contri[key] = value
                    else:
                        src_contri[key] = value 

        #Rank the node with highest src influence
        src_contri = dict(sorted(src_contri.items(), key=lambda item: item[1], reverse=True))
        #print("Source Contribution:\n%s" % (str(src_contri))) #Debug
        
        #Record stat for current generation
        fit_avg = round(sum(fit_list) / len(fit_list), 2)
        ga_stat.append([generation, fit_max, fit_avg, fit_set])

    #Generate the dataframe to show the result
    result_ga = pd.DataFrame(ga_stat, columns=['Generation', 'Best Fitness', 'Avg Fitness', 'Best Set'])
    result_ga.plot(x="Generation", y=["Best Fitness", "Avg Fitness"], kind="line")

    result_ga.to_csv(str(output_file) + '_GA.csv') #Generate the output
    plt.savefig(str(output_file) + '_GA.png') #Store the output graph

    end_time = time.time()
    total_time = end_time - start_time

    print ("GA end: " + str(output_file))
    return fit_set, fit_max, output_file, total_time

#########################################################

#Run GA Multiple Times
def run_GA(graph_file, inf, popu, chromo, max_gen, times, time_sir):
    if __name__ == "__main__":
        start_time = time.time()
        data = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(GeneAlgo, graph_file, inf, popu, chromo, max_gen, num, time_sir) for num in range(times)]
            for future in concurrent.futures.as_completed(futures):
                n_set, n_fit, name, t = future.result()
                data.append([name, n_set, n_fit, t])
        final_ga = pd.DataFrame(data, columns=['run', 'Fittest Set', 'Fitness', 'Time'])
        final_ga.to_csv(str(graph_file).strip(".txt") + '_GA.csv')

#Run SIR Multiple Times
def run_SIR(graph_file, init_SP, time_sir, inf, rec, times):
    if __name__ == "__main__":  
        data = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(sirModel, graph_file, init_SP, time_sir, inf, rec, str(num)) for num in range(times)]
            for future in concurrent.futures.as_completed(futures):
                run, infected = future.result()
                data.append([run, infected])
        final_sir = pd.DataFrame(data, columns=['run', 'Inf + Rec'])
        final_sir.to_csv(str(graph_file).strip(".txt") + '_' + str(inf) + '_SIR.csv')



#############################################
#Dataset file format: edge list text file (e.g. 0 1, edge from 0 to 1)
#Node should begin with 0 and node number inside the graph should be consecutive (i.e. 1,2,3,...)
s_fb = "facebook_combined.txt"
s_twitch = "musae_PTBR_edges.txt"
s_fb_food = "fb_pages_food_edges.txt"
s_hgdex = "hgd_example_graph.txt"
s_geekex = "GeeksforGeeks_example.txt"

src_ks_fb = [1912, 2543, 2347, 2266, 1985]
src_hgd_fb = [107, 2047, 1912, 1684, 348]
src_ga1_fb = [107, 754, 1919, 3631, 4025]
src_ga2_fb = [107, 698, 779, 2358, 3260]
src_ga3_fb = [107, 809, 1154, 1957, 3849]

src_ks_tw = [127, 1476, 290, 1297, 467]
src_hgd_tw = [127, 224, 208, 1188, 1421]
src_ga1_tw = [188, 467, 689, 1476, 1593]
src_ga2_tw = [67, 404, 486, 1374, 1734]
src_ga3_tw = [127, 290, 1297, 1455, 1543]

src_ks_pa = [265, 67, 340, 90, 56]
src_hgd_pa = [265, 518, 374, 555, 299]
src_ga1_pa = [31, 119, 180, 265, 340]
src_ga2_pa = [48, 246, 265, 441, 598]
src_ga3_pa = [119, 217, 317, 340, 508]

#SIR Simulation
#inf = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25]
#for i in inf:   
#    run_SIR(s_twitch, src_ga3_tw, 5, i, 1, 50)

graph = snap.LoadEdgeList(snap.TUNGraph, s_fb, 0, 1)

diff_nnode = set()
for src in src_ga3_fb:
    NI = graph.GetNI(src)
    nlist = [e for e in NI.GetOutEdges()]
    for n in nlist:
        diff_nnode.add(n)
print(len(list(diff_nnode)))

#Average shortest degree
#graph = snap.LoadEdgeList(snap.TUNGraph, s_fb, 0, 1)
#path = list(combinations(src_ks_tw, 2))
#dist = 0

#for c in path:
#    dist += graph.GetShortPath(c[0],c[1])

#print(dist/len(path))



#K-shell - Input: Dataset, no of spreader; Output: ks value dataframe(used in HGD), spreaders
#ks_df, ks_src = kshellDecomp(s_geekex, 3)


#HGD - Input: Dataset, ks value dataframe, no of spreader; Output: spreaders
#ks_df, ks_src = kshellDecomp(s_fb_food, 3)
#hgd_src = hgdAlgo(s_hgdex, ks_df, 3)
#print (str(hgd_src))


#GA - Input: Dataset, Infection Rate, Population, no of spreader, Generation, Output File Name, SIR Time (Evaluation)
#Output: Fittest node set, Fitness, Output File Name, Running Time
#src, fit, name, t = GeneAlgo(s_fb_food, 0.3, 12, 5, 3, "Demo", 3)
#print(str(src) + " " + str(fit))


#SIR Model - Input: Dataset, Source Node, Running Time (Iteration), Infection Rate (0 - 1), Recovery Rate (0 - 1), Output File Name
#Output: Output File Name, Infected Node (Infected + Recovered) at latest iteration (time)
#name, inf_num = sirModel(s_fb_food, hgd_src, 10, 0.5, 0.3, "Demo")


#Run GA Multiple times: Input: Dataset, Infection Rate, Population, no of spreader, Generation, GA Times, SIR Time (Evaluation)
#run_GA(s_fb, 0.5, 12, 5, 500, 10, 3)


#Run SIR Multiple times: Input: Dataset, Source Node, Running Time (Iteration in a SIR), Infection Rate, Recovery Rate, SIR Run





##############################################################














