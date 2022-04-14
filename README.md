# Identifying-Influential-Spreaders
Algorithms for finding influrencers with simulating model

There are 3 algorithms for finding influential spreaders
* K-Shell Decomposition [1]
* Heuristic Group Discovery (HGD) Algorithm [2]
* Genetic Algorithm [3]
- - - -
The **Suspectible-Infection-Recovery (SIR) Model [4]** is used to simulate the influrencing (spreading) process.
- - - - 
*Reference:*

[1]: “K-shell decomposition on Social Networks,” GeeksforGeeks, 01-Oct-2020. [Online]. Available: https://www.geeksforgeeks.org/k-shell-decomposition-on-social-networks/.

[2]: L. Jiang, X. Zhao, B. Ge, W. Xiao, and Y. Ruan, “An efficient algorithm for mining a set of influential spreaders in complex networks,” Physica A: Statistical Mechanics and its Applications, vol. 516, pp. 58–65, 2019.

[3]: K. C. Wong and K. Y. Szeto, “Multiple Sources Influence Maximization in Complex Networks with Genetic Algorithm,” Distributed Computing and Artificial Intelligence, 16th International Conference, pp. 226–234, 2019.

[4]: N. Antulov-Fantulin, A. Lančić, H. Štefančić, and M. Šikić, “FastSIR algorithm: A fast algorithm for the simulation of the epidemic spread in large networks by using the susceptible–infected–recovered compartment model,” Information Sciences, vol. 239, pp. 226–240, 2013.

- - - -
**Packages Used**
* snap: Access graphs
* pandas: Organize datas as tables.
* matplotlib: Visualize the result
* random: For creating random number to simulating infection process
* collections: Call the "Set" data structure
* concurrent: Run the repeated function with parallelism to reduce total running time
* time: Record the running time of the algorithms
