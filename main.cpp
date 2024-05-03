// Common CGP implementation (cgp.cpp)
// ========================================
// Cartesian genetic programming, introduced by Julian Miller and Peter Thomson
// in 2000, is a variant of genetic programming where the genotype is represented
// as a list of integers that are mapped to directed oriented graphs.
// In order to evaluate the fitness function the response for each training vector
// has to be calculated. This step involves the interpretation of a CGP genotype
// for each vector. In this implementation we are using *linear interpreter* that
// interprets CGP chromosome.
//
/*===============================================================================
  cgp.cpp: Symbolic regression in floating-point domain with CGP
  ===============================================================================
  Copyright (C) 2012 Brno University of Technology,
                     Faculty of Information Technology

  Author(s): Zdenek Vasicek <vasicek AT fit.vutbr.cz>
  =============================================================================*/

/// #define HAVE_POPCNT
/// #define DONOTEVALUATEUNUSEDNODES
#define MAX_POPSIZE 100
typedef double fitvaltype;
typedef double nodetype;


#include "cgp.h"
#include "symb.h"
#include "math.h"
#include "cgpROC.cpp"
#include <vector>

pchromosome population[MAX_POPSIZE];
fitvaltype fitvalues[MAX_POPSIZE];
/// default parameters
tparams params = {50000 /*generations*/, 5 /*pop.size*/, 5 /*mut.genes*/, 0, 0, 15 /*cols*/, 1 /*rows*/, 15 /*lback*/,
                  2, 1, 9 /*functions*/};

nodetype *nodeoutput; //array of node outputs used for CGP evaluation
int *isused[MAX_POPSIZE]; //array of marked used nodes for each individual
int usednodes[MAX_POPSIZE]; //number of used nodes for each individual
fitvaltype *data;  //training data
col_validvalues **col_values; //valid gene values for each column

// CGP mutation operator
// ----------------------------------------------------------------------------
// This function modifies `params.mutations` randomly selected genes of a given genotype. The modification is done in place.
//
// The implementation of the mutation operator has to ensure that the modifications are legal and lead to a valid phenotype.
// This is done using `col_values` array which contains valid values that can occur in each column of CGP array.

inline void cgp_mutate(chromosome p_chrom) {
    int rnd = rand() % params.mutations;
    int genes = rnd + 1;
    for (int j = 0; j < genes; j++) {
        int i = rand() % (params.geneoutidx + params.outputs);
        int col = (int) (i / params.genespercol);

        if (i < params.geneoutidx) {
            if ((i % params.nodeios) < params.nodeinputs) {  //mutate a gene that encodes connection
                if (col_values[col]->items > 1) {
                    do { rnd = col_values[col]->values[(rand() % (col_values[col]->items))]; }
                    while (rnd == p_chrom[i]);
                    p_chrom[i] = rnd;
                }
            } else {  //mutate a gene that encodes function
                if (params.nodefuncs > 1) {
                    do { rnd = rand() % params.nodefuncs; } while (rnd == p_chrom[i]);
                    p_chrom[i] = rnd;
                }
            }
        } else {  //mutate a primary output
            do { rnd = rand() % params.lastnodeidx; } while (rnd == p_chrom[i]);
            p_chrom[i] = rnd;
        }
    }
}

// CGP evaluation
// ----------------------------------------------------------------------------
// This function simulates the encoded candidate solution and calculates
// response for single input vector. The response for each CGP node is
// stored in array `nodeoutput`. The first `params.inputs` items contain
// input vector.
//
// The interpreter consists of a loop that calculates the response for each CGP
// node according to the genes of CGP chromosome. The execution of the encoded graph
// starts from the first node and continues according to the increasing node index.
// This scheme represents the most efficient
// implementation as it does not introduce any overhead due to function calling that
// have to manipulate with stack. In contrast with the recursive approach, all the
// output values are calculated in one pass. In order to improve the performance,
// a simple preprocessing step that marks the utilized nodes only can be introduced
// (by enabling `DONOTEVALUATEUNUSEDNODES` parameter). Only the marked nodes are
// subsequently evaluated.
//
inline void cgp_eval(chromosome p_chrom, int *isused) {
    int fce, i, j;
    nodetype in1, in2;
    nodetype *pnodeout = nodeoutput + params.inputs;

#ifdef DONOTEVALUATEUNUSEDNODES
    isused += params.inputs;
#endif

    /// Evaluate the response of each node
    for (i = 0; i < params.cols; i++)
        for (j = 0; j < params.rows; j++) {
#ifdef DONOTEVALUATEUNUSEDNODES
            if (!*isused++) {  //This node is not used, skip it
                p_chrom += 3;
                pnodeout++;
                continue;
            }
#endif

            in1 = nodeoutput[*p_chrom++];
            in2 = nodeoutput[*p_chrom++];
            fce = *p_chrom++;
            switch (fce) {
                case 0:
                    *pnodeout++ = 0.25;
                    break;         //const
                case 1:
                    *pnodeout++ = 0.50;
                    break;         //const
                case 2:
                    *pnodeout++ = 1.00;
                    break;         //const
                case 3:
                    *pnodeout++ = in1;
                    break;          //in1
                case 4:
                    *pnodeout++ = in1 + in2;
                    break;    //in1 + in2
                case 5:
                    *pnodeout++ = in1 - in2;
                    break;    //in1 - in2
                case 6:
                    *pnodeout++ = in1 * in2;
                    break;    //in1 * in2
                case 7:
                    *pnodeout++ = (in2 == 0) ? 1e10 : in1 / in2;
                    break;     //in1 / in2
                case 8:
                    *pnodeout++ = sin(in1);;
                    break;   //sin(in1)
                default:
                    abort();
            }
        }
}

// Calculate fitness values
// ---------------------------------------------
// This function determines the fitness value of each candidate solution except
// parental solution `parentidx`. The calculated fitness values are stored
// in array denoted as `fitvalues`.
//
// In order to calculate the fitness value of each candidate solution the
// response for each training data has to be evaluated. In this implementation
// the fitness value corresponds with the Hamming distance between the calculated and
// desired response.
//
// The fitness function calls `cgp_eval` procedure that, given the CGP calculated outputs
// data and a candidate solution, evaluates the response for a single input vector
// (i.e. single fitness case stored in array denoted as `ptraindata`). Then, according to
// the information about output connections stored in chromosome, the fitness value of a
// genotype for the utilized input vector is calculated. These steps are repeated until last
// training vector is evaluated.
//
inline void calc_fitness(int parentidx) {
    chromosome p_chrom;
    fitvaltype fit, vysl;

    for (int i = 0; i < params.popsize; i++) {
        if (i == parentidx) {
            continue;
        }

#ifdef DONOTEVALUATEUNUSEDNODES
        usednodes[i] = used_nodes((chromosome) population[i], isused[i]);
#endif

        fitvalues[i] = 0;
    }

    ///determine and check response of each candidate solution
    for (int i = 0; i < params.popsize; i++) {
        nodetype *ptraindata = data;
        //printf("Candidate solution: %-8d\n", i);
        if (i == parentidx) continue;

        float tp = 0;
        float tn = 0;
        float fp = 0;
        float fn = 0;
        float accuracy;

        std::vector<fitvaltype> desired, obtained;


        for (int l = 0; l < params.trainingvectors; l++) {

            ///copy the first part of a training vector to the primary inputs
            memcpy(nodeoutput, ptraindata, params.inputs * sizeof(nodetype));
            ptraindata += params.inputs;

            cgp_eval((chromosome) population[i], isused[i]);

            ///compute sum of absolute differences between the calculated and desired output values
            fit = 0;

            p_chrom = (chromosome) population[i] + params.geneoutidx;


            for (int k = 0; k < params.outputs; k++) {

                desired.push_back(*(ptraindata + k));
                obtained.push_back(nodeoutput[*p_chrom++]);

            }

            ///next training vector
            ptraindata += params.outputs;
        }

        /// specify threshold for classification
        double threshold = get_threshold(desired, obtained);

        int class_obtained;
        for (int i = 0; i < obtained.size(); i++) {
            if (obtained[i] < threshold) {
                class_obtained = 0;
            } else {
                class_obtained = 1;
            }

            if (class_obtained == 1) {
                if (desired[i] == 1) {
                    tp += 1;
                } else {
                    fp += 1;
                }
            } else {
                if (desired[i] == 1) {
                    fn += 1;
                } else {
                    tn += 1;
                }
            }
        }

        /// calculate accuracy of candidate solution
        if (tp + fp == 0) {
            accuracy = 0;
        } else {
            accuracy = tp / (tp + fp);
        }
        fitvalues[i] = accuracy;
    }
}

// Number of utilized nodes
// ----------------------------------------------------------------------------
// This function calculates the number of nodes that contribute to the
// resulting phenotype that is encoded by a given genotype.
//
int used_nodes(chromosome p_chrom, int *isused) {
    int in, fce, idx, used = 0;
    int *pchrom;

    memset(isused, 0, params.lastnodeidx * sizeof(int));

    //mark nodes connected to the primary outputs
    pchrom = p_chrom + params.geneoutidx;
    for (int i = 0; i < params.outputs; i++)
        isused[*pchrom++] = 1;

    //go throught the cgp array
    pchrom = p_chrom + params.geneoutidx - 1;
    idx = params.lastnodeidx - 1;
    for (int i = params.cols; i > 0; i--) {
        for (int j = params.rows; j > 0; j--, idx--) {
            fce = *pchrom--; //fce
            if (isused[idx] == 1) {
                // the current node is marked, mark also the connected nodes
                in = *pchrom--; // in2
                if (fce > 3 && fce < 8)    // 2-input functions
                    isused[in] = 1;
                in = *pchrom--; // in1
                if (fce > 2)
                    isused[in] = 1;

                used++;
            } else {
                // the current node is not market, skip it
                pchrom -= params.nodeinputs;
            }
        }
    }

    return used;
}

std::string *get_eq(chromosome p_chrom, int nodeidx) {
    using namespace std;
    string *s = new string();
    char buf[50];
    if (nodeidx < params.inputs) {
        sprintf(buf, "in_%d", nodeidx);
        *s += (nodeidx == 0) ? "x" : (nodeidx == 1) ? "y" : (nodeidx == 2) ? "z" : buf;
    } else {
        string *in1 = get_eq(p_chrom, *(p_chrom + (params.nodeios) * (nodeidx - params.inputs)));
        string *in2 = get_eq(p_chrom, *(p_chrom + (params.nodeios) * (nodeidx - params.inputs) + 1));
        int fce = *(p_chrom + (params.nodeios) * (nodeidx - params.inputs) + params.nodeios - 1);

        switch (fce) {
            case 0:
                *s = "0.25";
                break;
            case 1:
                *s = "0.50";
                break;
            case 2:
                *s = "1.00";
                break;
            case 3:
                *s = *in1;
                break;                               //in1
            case 4:
                *s = "(" + *in1 + " + " + *in2 + ")";
                break;    //in1 + in2
            case 5:
                *s = "(" + *in1 + " - " + *in2 + ")";
                break;    //in1 - in2
            case 6:
                *s = "(" + *in1 + " * " + *in2 + ")";
                break;    //in1 * in2
            case 7:
                *s = "(" + *in1 + " / " + *in2 + ")";
                break;    //in1 / in2
            case 8:
                *s = "sin(" + *in1 + ")";
                break;                //sin(in1)
            default:
                abort();
        }
        delete in1, in2;
    }

    return s;
}

void print_eq(FILE *fout, chromosome p_chrom) {
    using namespace std;

    for (int i = 0; i < params.outputs; i++) {
        string *s = get_eq(p_chrom, *(p_chrom + params.geneoutidx + i));
        fprintf(fout, "%s", s->c_str());
        delete s;
    }
}

// Main application
// -------------------------------------
// Symbolic regression problem using Cartesian Genetic Programming.
//
// Usage:
//
//     cgp data.txt [-c COLUMNS] [-r ROWS] [-l LBACK] \
//     [-g GENERATIONS] [-p POPSIZE] [-m MUTATEDGENES] [-e max error]
//
int main(int argc, char *argv[]) {
    using namespace std;
    int blk, bestblk, data_items, parentidx, fittest_idx;
    unsigned long int generation;
    fitvaltype bestfitval;

    params.accfitval = 0;
    strcpy(params.datafname, "banknotes.txt");
    parse_options(argc, argv);
    ///load training data
    if ((data_items = parsefile(params.datafname, NULL, NULL, NULL)) < 1) {
        printf("Invalid data\n");
        return 0;
    }
    data = new nodetype[data_items];
    parsefile(params.datafname, data, &params.inputs, &params.outputs);
    params.trainingvectors = data_items / (params.inputs + params.outputs); //Spocitani poctu pruchodu pro ohodnoceni
    printf("Training data:\n   file:%s, vectors:%d, inputs:%d, outputs:%d\n", params.datafname, params.trainingvectors,
           params.inputs, params.outputs);

    init_paramsandluts();

    ///memory allocation
    for (int i = 0; i < params.popsize; i++) {
        population[i] = new chromosome[params.geneoutidx + params.outputs];
        isused[i] = new int[params.genes];
    }
    nodeoutput = new nodetype[params.genes];

    srand((unsigned) time(NULL));

    // Initial population and its evaluation
    // -----------------------------------------------------------------------
    // The initial population which consists of `params.popsize` individuals
    // is generated randomly.

    for (int i = 0; i < params.popsize; i++) {
        chromosome p_chrom = (chromosome) population[i];
        for (int j = 0; j < params.cols * params.rows; j++) {
            int col = (int) (j / params.rows);
            for (int k = 0; k < params.nodeinputs; k++)  // node inputs
                *p_chrom++ = col_values[col]->values[(rand() % (col_values[col]->items))];
            *p_chrom++ = rand() % params.nodefuncs; // node function
        }
        for (int j = 0; j < params.outputs; j++) // primary outputs
            *p_chrom++ = rand() % params.lastnodeidx;
    }

    // Evaluate the initial population and find the fittest candidate solution
    // that becomes parent
    calc_fitness(-1);
    bestfitval = fitvalues[0];
    fittest_idx = 0;
    for (int i = 1; i < params.popsize; i++)
        if (fitvalues[i] < bestfitval) {
            bestfitval = fitvalues[i];
            fittest_idx = i;
        }

    printf("CGP parameters:\n   l-back=%d, rows=%d, cols=%d, functions=%d\n", params.lback, params.rows, params.cols,
           params.nodefuncs);
    printf("   popsize=%d, mutgenes=%d, generations=%d, acceptable error=%f\n", params.popsize, params.mutations,
           params.maxgenerations, params.accfitval);
#ifdef DONOTEVALUATEUNUSEDNODES
    printf("   evalunused=false\n");
#endif
    printf("Evolutionary run:\n   initial fitness=%f\n", bestfitval);
    double time = cpuTime();

    // Evolutionary loop
    // -----------------------------------------------------------------------
    for (int generation = 0; generation < params.maxgenerations; generation++) {
        //printf("Generation: %-8d\n",generation);
        // ### Step 1 ###
        // Generate offsprings of the fittest individual
        for (int i = 0; i < params.popsize; i++) {
            if (fittest_idx == i) continue;
            cgp_mutate(copy_chromozome(population[fittest_idx], population[i]));
        }

        // ### Step 2 ###
        // Evaluate population and calculate fitness values
        calc_fitness(fittest_idx);
        // ### Step 3 ###
        // Check if there is an offspring that should replace parental solution
        int newparidx = fittest_idx;
        for (int i = 0; i < params.popsize; i++) {
            fitvaltype fit = fitvalues[i];
            if ((i == fittest_idx) || (bestfitval < fit)) continue;

            // current individual is at least of the same quality as its parent
            if (fit > params.accfitval) {  //current individual does not
                if (fit < bestfitval) {
                    printf("Generation: %-8d\tFitness: %f\n", generation, fit);
                    //print_chrom(stdout, (chromosome)population[i]);
                    //printf("Equation: "); print_eq(stdout, (chromosome)population[i]); printf("\n");
                }
                bestfitval = fit;
                bestblk = params.nodes;
                newparidx = i;
            } else {  //second optimization criterion - the number of utilized nodes
#ifdef DONOTEVALUATEUNUSEDNODES
                blk = usednodes[i];
#else
                blk = used_nodes((chromosome) population[i], isused[i]);
#endif
                if (blk <= bestblk) {
                    if ((blk < bestblk) || (fit < bestfitval)) {
                        //printf("Generation: %-8d\tFitness: %f\tNodes: %d\n",generation, fit, blk);
                        //printf("Equation: "); print_eq(stdout, (chromosome)population[i]); printf("\n");
                    }

                    bestfitval = fit;
                    bestblk = blk;
                    newparidx = i;
                }
            }

        }
        fittest_idx = newparidx;
    }

    // End of evolution
    // -----------------------------------------------------------------------
    time = cpuTime() - time;
    printf("Best fitness: %f\n", bestfitval);
    printf("Best individual: ");
    print_chrom(stdout, (chromosome) population[fittest_idx]);
    printf("Duration: %f Evaluations per sec: %f\n", time, params.maxgenerations * (params.popsize - 1) / time);

    if (bestfitval < 0.5) {
        //save the solution
        FILE *chrfil = fopen("result.chr", "wb");
        print_chrom(chrfil, (chromosome) population[fittest_idx]);
        fclose(chrfil);
    }

    return 0;
}
