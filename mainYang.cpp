//
//  main.cpp
//  Community Structure Code (Yang et al approach)
//
//  Created by Theodore Faust on 10/7/23.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <stdexcept>
#include <vector>
#include "MathMatrix.h"
#include "MathVector.h"
#include <algorithm>
#ifdef _MSC_VER
#include "iso646.h"      // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#endif
std::vector<MathMatrix> gibbsStatistics(MathMatrix G, int numGroups, std::vector<MathMatrix> A);
long double logobjFunc(MathMatrix G, long double temperature, int numGroups, std::vector<MathMatrix> A);
MathVector objFuncProbs(MathMatrix G, int node, int layer, int numGroups, long double temperature, std::vector<MathMatrix> A);
int NodeSample(MathMatrix G, int node, int layer, int numGroups, long double temperature, std::vector<MathMatrix> A);
long double lbeta(long double x, long double y);
int main(int argc, const char * argv[])
{
    std::random_device rd3;
    std::mt19937 gen3(rd3());
    std::uniform_real_distribution<> dis3(0.0, 1.0);
    
    int numGroups = 2;
    size_t numNodes = 100;
    size_t sizeA = numNodes;
    std::vector<int> group1Sizes = {50,60,70,80,90};
    size_t numRuns = group1Sizes.size();
    std::vector<int> offsets = {0,-5,-10,-5,0};
    size_t layersA = offsets.size();
    size_t numLayers = layersA;
    int runs = 500;
    std::vector<double> sigmaList = {0.25};
    double perturbProb = 3e-5;
    int numPerturbSteps = 50;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0,numLayers-1);
    std::random_device rd2;
    std::mt19937 gen2(rd2());
    std::uniform_int_distribution<> distrib2(0,numGroups-1);
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::uniform_int_distribution<> distrib4(0,numNodes-1);
    long numMCSteps = 1e3;
    bool runSingleLayer = false;
    bool runMultiLayer = true;
    bool noGpSizeBias = true;
    std::vector<double> alphaList = {1};
    std::vector<double> betaList = {1};
    
    if (runMultiLayer)
    {
        for (size_t sigmaIndex = 0; sigmaIndex < sigmaList.size(); ++ sigmaIndex)
        {
            std::vector<std::string> stepTypes = {"yang"};
            for (size_t stepType = 0; stepType < 1; stepType++)
            {
                for (size_t group1Index = 0; group1Index < numRuns; ++group1Index)
                {
                    std::vector<MathVector> numCorrect(numLayers,MathVector(runs));
                    std::vector<MathMatrix> initialGroups(numLayers,MathMatrix(numNodes,runs));
                    std::vector<MathMatrix> finalGroups(numLayers,MathMatrix(numNodes,runs));
                    std::vector<MathMatrix> Gcount(numGroups,MathMatrix(numNodes,numLayers));
                    MathMatrix sbmProbs = {{sigmaList[sigmaIndex],0.1},{0.1,sigmaList[sigmaIndex]}};
                    std::vector<MathMatrix> A(layersA,MathMatrix(sizeA,sizeA));
                    size_t group1Size = group1Sizes[group1Index];
                    size_t group0Size = numNodes-group1Size;
                    int groupI = 0;
                    int groupJ = 0;
                    MathVector iterationNumbers = {20,10,10,10,10,10,10,5,5,5};
                    MathVector temperatures = {1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
                    size_t gibbsSteps = iterationNumbers.size();
                    for (size_t l = 0; l < layersA; ++l)
                    {
                        for (size_t i = 0; i < sizeA; ++i)
                        {
                            if (i < (group1Size + offsets[l]))
                            {
                                groupI = 1;
                            }
                            else
                            {
                                groupI = 0;
                            }
                            for (size_t j = 0; j < i; ++j)
                            {
                                if (j < (group1Size + offsets[l]))
                                {
                                    groupJ = 1;
                                }
                                else
                                {
                                    groupJ = 0;
                                }
                                double randomNum = dis3(gen3);
                                if (randomNum < sbmProbs(groupI,groupJ))
                                {
                                    A[l](i,j) = 1.0;
                                    A[l](j,i) = 1.0;
                                }
                                else
                                {
                                    A[l](i,j) = 0.0;
                                    A[l](j,i) = 0.0;
                                }
                            }
                        }
                    }
                    for (size_t l = 0; l < layersA; ++l)
                    {
                        std::ofstream myfile3("layer" + std::to_string(l) + "group1size" + std::to_string(group1Size) + "sigma" + std::to_string(sigmaList[sigmaIndex]) + ".txt",std::ios::trunc);
                        myfile3 << A[l];
                        myfile3.close();
                    }
                    for (size_t alphaIndex = 0; alphaIndex < alphaList.size(); ++alphaIndex)
                    {
                        for (size_t betaIndex = 0; betaIndex < betaList.size(); ++betaIndex)
                        {
                            double alpha = alphaList[alphaIndex];
                            double beta = betaList[betaIndex];
                            for (size_t run = 0; run < runs; ++run)
                            {
                                bool stoppingCond = false;
                                int stepCounter = 0;
                                MathMatrix G = MathMatrix(numNodes,numLayers);
                                for (size_t i = 0; i < numNodes; ++i)
                                {
                                    for (size_t l = 0; l < numLayers; ++l)
                                    {
                                        G(i,l) = distrib2(gen2);
                                    }
                                }
                                for (size_t i = 0; i < numNodes; ++i)
                                {
                                    for (size_t l = 0; l < numLayers; ++l)
                                    {
                                        initialGroups[l](i,run) = G(i,l);
                                    }
                                }
                                MathMatrix Gtemp = G;
                                for (size_t step = 0; step < gibbsSteps; ++step)
                                {
                                    for (size_t innerStep = 0; innerStep < iterationNumbers(step); ++innerStep)
                                    {
                                        for (size_t i = 0; i < numNodes; ++i)
                                        {
                                            for (size_t l = 0; l < numLayers; ++l)
                                            {
                                                int newGroup = NodeSample(Gtemp, i, l, numGroups, temperatures(step), A);
                                                Gtemp(i,l) = newGroup;
                                            }
                                        }
                                    }
                                }
                                G = Gtemp;
                                std::cout << "run " << run << " complete" << "\n";
                                for (size_t l = 0; l < numLayers; ++l)
                                {
                                    int numCorrectRun = 0;
                                    for (size_t q = 0; q < numNodes; ++q)
                                    {
                                        finalGroups[l](q,run) = G(q,l);
                                        if (q < group1Size + offsets[l])
                                        {
                                            if (fabs(G(q,l) - 1.0) < 0.1)
                                            {
                                                numCorrectRun++;
                                            }
                                        }
                                        else
                                        {
                                            if (fabs(G(q,l)) < 0.1)
                                            {
                                                numCorrectRun++;
                                            }
                                        }
                                    }
                                    
                                    if (numNodes - numCorrectRun > numCorrectRun)
                                    {
                                        numCorrectRun = numNodes - numCorrectRun;
                                    }
                                    numCorrect[l](run) = numCorrectRun;
                                }
                                for (size_t l = 0; l < numLayers; ++l)
                                {
                                    std::ofstream myfile3("layer" + std::to_string(l) + "numCorrect" + std::to_string(group1Size) + stepTypes[stepType] + "sigma" + std::to_string(sigmaList[sigmaIndex]) + "new.txt",std::ios::trunc);
                                    myfile3 << numCorrect[l];
                                    myfile3.close();
                                    std::ofstream myfile4("layer" + std::to_string(l) + "initialGroups" + std::to_string(group1Size) + stepTypes[stepType] + "sigma" + std::to_string(sigmaList[sigmaIndex]) +  "new.txt",std::ios::trunc);
                                    myfile4 << initialGroups[l];
                                    myfile4.close();
                                    std::ofstream myfile5("layer" + std::to_string(l) + "finalGroups" + std::to_string(group1Size) + stepTypes[stepType] + "sigma" + std::to_string(sigmaList[sigmaIndex]) +  "new.txt",std::ios::trunc);
                                    myfile5 << finalGroups[l];
                                    myfile5.close();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

long double mcIntegrand(MathVector kronDel1, MathVector kronDel2, double p1, double p2)
{
    if (kronDel1.size() != kronDel2.size())
    {
        throw std::out_of_range("Sizes of vectors unequal");
    }
    else if (p1 > 1.0 || p1 < 0.0 || p2 > 1.0 || p2 < 0.0)
    {
        throw std::out_of_range("Probability not between 0 and 1");
    }
    else
    {
        size_t kdSize = kronDel1.size();
        long double out = 1.0;
        for (size_t i = 0; i < kdSize; ++i)
        {
            out *= (kronDel1(i)*p1 + (1.0-p1)*(kronDel2(i)*p2 + (1.0-kronDel2(i))*(1.0-p2)));
        }
        return out;
    }
}

int NodeSample(MathMatrix G, int node, int layer, int numGroups, long double temperature, std::vector<MathMatrix> A)
{
    MathVector probList = objFuncProbs(G, node, layer, numGroups, temperature, A);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double tempRand = dis(gen);
    int tempGroup = 0;
    while (tempRand > probList(tempGroup))
    {
        tempGroup++;
    }
    return tempGroup;
}

MathVector objFuncProbs(MathMatrix G, int node, int layer, int numGroups, long double temperature, std::vector<MathMatrix> A)
{
    MathMatrix Gtemp = G;
    MathVector temp = MathVector(numGroups);
    MathVector exptemp = MathVector(numGroups);
    MathVector out = MathVector(numGroups);
    for (size_t g = 0; g < numGroups; ++g)
    {
        Gtemp(node,layer) = g;
        long double prob = logobjFunc(Gtemp, temperature, numGroups, A);
        temp(g) = prob;
    }
    double maxLogProb = -INFINITY;
    double secondMaxLogProb = -INFINITY;
    size_t maxLogProbIndex = 0;
    for (size_t g = 0; g < numGroups; ++g)
    {
        if (temp(g) > maxLogProb)
        {
            maxLogProb = temp(g);
            maxLogProbIndex = g;
        }
        else if (temp(g) > secondMaxLogProb)
        {
            secondMaxLogProb = temp(g);
        }
    }
    if ((maxLogProb - secondMaxLogProb) > 35)
    {
        exptemp(maxLogProbIndex) = 1;
    }
    else
    {
        temp -= MathVector(numGroups,maxLogProb);
        for (size_t g = 0; g < numGroups; ++g)
        {
            exptemp(g) = std::exp(temp(g));
        }
    }
    exptemp = exptemp/exptemp.norm1();
    for (size_t g = 0; g < numGroups; ++g)
    {
        for (size_t i = g; i < numGroups; ++i)
        {
            out(i) += exptemp(g);
        }
    }
    return out;
}

long double logobjFunc(MathMatrix G, long double temperature, int numGroups, std::vector<MathMatrix> A)
{
    std::vector<MathMatrix> stats = gibbsStatistics(G, numGroups, A);
    double logOut = 0.0;
    for (size_t k = 0; k < numGroups; ++k)
    {
        logOut += std::lgamma(stats[0](k,1) + 1);
    }
    for (size_t k = 0; k < numGroups; ++k)
    {
        for (size_t l = 0; l < numGroups; ++l)
        {
            logOut += std::lgamma(stats[3](k,l) + 1);
        }
        logOut -= std::lgamma(stats[4](k,1) + numGroups);
    }
    for (size_t k = 0; k < numGroups; ++k)
    {
        for (size_t l = k+1; l < numGroups; ++l)
        {
            logOut += lbeta(stats[2](k,l) + 1,stats[1](k,l) - stats[2](k,l) + 1);
        }
    }
    for (size_t k = 0; k < numGroups; ++k)
    {
        logOut += lbeta((stats[2](k,k))/2 + 1, (stats[1](k,k)-stats[2](k,k))/2 + 1);
    }
    return logOut/temperature;
}

std::vector<MathMatrix> gibbsStatistics(MathMatrix G, int numGroups, std::vector<MathMatrix> A)
{
    std::vector<MathMatrix> out(5,MathMatrix(numGroups,numGroups));
    size_t numNodes = G.getRowSize();
    size_t numLayers = G.getColSize();
    for (size_t i = 0; i < numNodes; ++i)
    {
        out[0](static_cast<int>(G(i,1)),1) = out[0](static_cast<int>(G(i,1)),1) + 1;
    }
    for (size_t l = 0; l < numLayers; ++l)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            for (size_t j = 0; j < numNodes; ++j)
            {
                out[1](static_cast<int>(G(i,l)),static_cast<int>(G(j,l))) = out[1](static_cast<int>(G(i,l)),static_cast<int>(G(j,l))) + 1;
                out[1](static_cast<int>(G(j,l)),static_cast<int>(G(i,l))) = out[1](static_cast<int>(G(j,l)),static_cast<int>(G(i,l))) + 1;
                out[2](static_cast<int>(G(i,l)),static_cast<int>(G(j,l))) = out[2](static_cast<int>(G(i,l)),static_cast<int>(G(j,l))) + A[l](i,j);
                out[2](static_cast<int>(G(j,l)),static_cast<int>(G(i,l))) = out[2](static_cast<int>(G(j,l)),static_cast<int>(G(i,l))) + A[l](i,j);
            }
        }
    }
    for (size_t l = 1; l < numLayers; ++l)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            out[3](static_cast<int>(G(i,l-1)),static_cast<int>(G(i,l))) = out[3](static_cast<int>(G(i,l-1)),static_cast<int>(G(i,l))) + 1;
            out[4](static_cast<int>(G(i,l)),1) = out[4](static_cast<int>(G(i,l)),1) + 1;
        }
    }
    out[4] = out[4] + out[0];
    return out;
}

long double lbeta(long double x, long double y)
{
    return (std::lgamma(x)+std::lgamma(y))-(std::lgamma(x+y));
}
