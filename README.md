# Fitness Dependence Optimizer

The Fitness Dependence Optimizer (FDO) is a recently developed swarm intelligence algorithm used to solve optimization problems. It's inspired by the collective decision-making process during the bee swarming reproductive phase. Here's a breakdown of its key features:

Inspiration: Unlike other bee-inspired algorithms like Honey Bee Algorithm or Artificial Bee Colony Optimization, FDO draws inspiration from the actual swarming behavior of bees rather than their food foraging process.

PSO-based: Although not directly connected to Particle Swarm Optimization (PSO), FDO shares similarities in its update mechanism for search agents' positions based on their velocity (called "pace" in FDO).

Key differentiator: What sets FDO apart is its reliance on fitness-dependent weights for each search agent. These weights dynamically adjust based on the agent's fitness value, influencing the balance between exploration (searching new areas) and exploitation (refining solutions in promising areas).

## Author of the Code Program
Jaza Mahmood Abdullah  jazamahmood@gmail.com

## Usage

```java
public static void main(String []args){
        FDOSingleObjective fdo = new FDOSingleObjective("FX1");
        fdo.setWeightFactor(0.0);
        fdo.setShouldPrintGlobalBestForEachTurn(true);
        fdo.setAgentRemainInsideBoundery(true);
        fdo.setNumberOfAgents(30);
        fdo.setMaxIteration(500);
        fdo.setTurns(1);
        fdo.runFDO();
        fdo.printMean();
        fdo.printStandardDeviation();
    }

```
