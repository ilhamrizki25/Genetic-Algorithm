
import static java.lang.Math.pow;
import java.util.Random;
import java.text.DecimalFormat;

/**
 * Genetic Algorithm By ilham Rizki Julianto - 1301170293 - IF-41-06
 */

class Genetic_Algorithm {

  static int populationSize = 50;
  static int maxGeneration = 100;
  static int chromosomeSize = 6; // Dont change this, since the chromosome design not dynamic
  static double probabilityCrossover = 0.80; // Range 0.1 to 1.0
  static double probabilityMutation = 0.07; // Range 0.1 to 1.0

  static Random rand = new Random();

  // Initialize population
  static double[][] initPopulation(int popSize, int chromosomeSize) {
    // Array for generated population stored
    double[][] generatedPopulation = new double[popSize][chromosomeSize];

    // Generate population
    for (int i = 0; i < popSize; i++) {
      for (int j = 0; j < chromosomeSize; j++) {
        generatedPopulation[i][j] = rand.nextDouble(); // Random real number
      }
    }

    return generatedPopulation;
  }

  // Print the population
  static void printPopulation(double[][] population, int chromosomeSize, double[] fitness) {
    for (int i = 0; i < population.length; i++) {
      System.out.print(i + ". ");
      for (int j = 0; j < chromosomeSize; j++) {
        System.out.printf("%f ", population[i][j]);
      }
      System.out.printf(" Fitness = %f", fitness[i]);
      System.out.println("");
    }
  }

  // Decode chromosome to define the value of x1 and x2
  static double[] decodeReal(double[] Chromosome) {

    // Boundaries of x1
    int max1 = 3;
    int min1 = -3;

    // Boundaries of x2
    int max2 = 2;
    int min2 = -2;

    // Formula to convert genotype to phenotype x1 and x2
    double x1 = min1 + ((max1 - min1) / 3) * (Chromosome[0] + Chromosome[1] + Chromosome[2]);
    double x2 = min2 + ((max2 - min2) / 3) * (Chromosome[3] + Chromosome[4] + Chromosome[5]);

    // Result of x1 and x2 saved to an array
    double[] result = new double[2];
    result[0] = x1;
    result[1] = x2;

    return result;
  }

  // Calculate objective value
  static double objective(double[] x1x2Value) {

    double x1 = x1x2Value[0];
    double x2 = x1x2Value[1];

    // Expression given in the paper
    double firstTerm = (4 - (2.1 * Math.pow(x1, 2)) + (Math.pow(x1, 4) / 3)) * Math.pow(x1, 2);
    double secondTerm = x1 * x2;
    double thirdTerm = (-4 + (4 * Math.pow(x2, 2))) * Math.pow(x2, 2);

    double result = firstTerm + secondTerm + thirdTerm;

    return result;
  }

  // Calculate ftiness value
  static double fitness(double objective) {
    return -objective;
  }

  // Calculate all fitness value for each chromosome in one generation
  static double[] allFitness(double[][] population) {
    double[] allFitness = new double[population.length];

    // Calculate all fitness value for each chromosome
    for (int i = 0; i < population.length; i++) {
      allFitness[i] = fitness(objective(decodeReal(population[i])));
    }

    return allFitness;
  }

  // Parent selection with roulette wheel selection method
  static double[][] rouletteWheel(double[][] population, double[] fitness) {
    double[] proportion = fitness.clone();
    double totalFitness = 0;
    double[][] parent = new double[10][];

    // Count the total of fitness point in one generation
    for (double fit : fitness) {
      totalFitness += fit;
    }

    // Calculate the value of proportion for each chromosome
    for (int i = 0; i < proportion.length; i++) {
      proportion[i] = fitness[i] / totalFitness;
    }

    // Pick 10 chromosome using roullete wheel method
    for (int i = 0; i < parent.length; i++) {
      double random = rand.nextDouble();
      double total = 0;
      int j = 0;
      while (total < random) {
        total += proportion[j];
        j++;
      }
      parent[i] = population[j - 1];
    }

    return parent;
  }

  // Parent selection with tournament selection method
  static double[][] tournament(double[][] population, double[] fitness) {
    double[][] temp = new double[4][];
    double[][] parent = new double[2][];
    double[] tempFitness = fitness.clone();

    // Select random candidate chromosome from population
    for (int k = 0; k < parent.length; k++) {
      for (int i = 0; i < temp.length; i++) {
        int random = rand.nextInt(population.length);
        temp[i] = population[random].clone();
        tempFitness[i] = fitness[random];
      }

      double bestFitness = tempFitness[0];
      int index = 0;

      // Select chromosome with the best fitness value from the candidate chromosome
      for (int j = 0; j < temp.length; j++) {
        if (bestFitness < tempFitness[j]) {
          index = j;
        }
      }

      // Selected candidate stored here
      parent[k] = temp[index].clone();
    }

    return parent;
  }

  // Crossover between child and parent to get the new generation of choromosome
  static double[][] crossover(double[][] parent, double probabilityCrossover) {
    double random = rand.nextDouble();

    if (random < probabilityCrossover) {
      double[][] child = parent.clone();

      // Crossover
      for (int j = 0; j < parent.length; j++) {
        int crossPoint = rand.nextInt(6);
        for (int i = crossPoint; i < 6; i++) {
          child[j][i] = parent[j + 1][i];
          child[j + 1][i] = parent[j][i];
        }
        j++;
      }
      return child;
    }

    return parent;
  }

  // Chromosome Mutation
  static double[] mutation(double[] child, double probabilityMutation) {
    for (int i = 0; i < 6; i++) {
      double random = rand.nextDouble();
      if (random < probabilityMutation) {
        child[i] = rand.nextDouble();
      }
    }

    return child;
  }

  // Survivor selection using general replacement
  static double[][] generalReplacement(double[][] oldPopulation, double[][] newPopulation, double[] fitness) {
    double maxFitness = fitness[0];
    int idx = 0;

    // Find the best fitness value
    for (int i = 0; i < fitness.length; i++) {
      if (maxFitness <= fitness[i]) {
        maxFitness = fitness[i];
        idx = i;
      }
    }

    // chromosome with best fitness value from old population copied 8 times to the new population with different placement
    for (int j = 0; j < 8; j++) {
      int random = rand.nextInt(fitness.length);
      newPopulation[random] = oldPopulation[idx].clone();
    }

    return newPopulation;
  }

  // Main Method
  public static void main(String[] args) {
    double[][] population = initPopulation(populationSize, chromosomeSize);
    double[] fitness;

    // printPopulation(population, chromosomeSize, fitness);

    for (int i = 0; i < maxGeneration; i++) {
      System.out.print("Generasi Ke-");
      System.out.println(i + 1);

      fitness = allFitness(population);
      // printPopulation(population, chromosomeSize, fitness);

      double[][] newPopulation = population.clone();
      int j = 0;
      while (j < populationSize) {
        double[][] parent = tournament(population, fitness);
        double[][] child = crossover(parent, probabilityCrossover);
        for (int k = 0; k < child.length; k++) {
          newPopulation[j++] = mutation(child[k], probabilityMutation).clone();
        }
      }

      // fitness = allFitness(newPopulation);
      // printPopulation(newPopulation, chromosomeSize, fitness);
      // fitness = allFitness(population);
      population = generalReplacement(population, newPopulation, fitness);
      fitness = allFitness(population);
      // printPopulation(population, chromosomeSize, fitness);

      double bestFitness = fitness[0];
      int idx = 0;
      for (int l = 0; l < fitness.length; l++) {
        if (bestFitness <= fitness[l]) {
          idx = l;
          bestFitness = fitness[l];
        }
      }

      double[] x1x2Value = decodeReal(population[idx]);

      // Outputing the result
      System.out.println("Found at Iterative  = " + idx);
      System.out.println("Value of x1         = " + x1x2Value[0]);
      System.out.println("Value of x2         = " + x1x2Value[1]);
      System.out.println("Fitness Value       = " + bestFitness);
      System.out.println("Objective Value     = " + objective(x1x2Value));
      System.out.println("-----------------------------------------------------");
    }
  }
}