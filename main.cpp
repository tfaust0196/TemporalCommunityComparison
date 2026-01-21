//
//  main.cpp
//  Community Structure Code
//
//  Created by Theodore Faust on 10/7/23.
//

#include <iostream>
#include <fstream>
#include <random>
#include <stdexcept>
#include <vector>
#include "MathMatrix.h"
#include "MathVector.h"
#include <algorithm>
#include <complex>
#include <chrono>
#include <cmath>
using namespace std::chrono;
using namespace std::complex_literals;
#ifdef _MSC_VER
#include "iso646.h"      // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#endif
long double logmcIntegrand(MathVector kronDel1, MathVector kronDel2, double p1, MathVector p2);
long double sumlog(long double x, long double y);
long double logbinom(double n, double k);
long double logCorrectionBazzi(MathMatrix G, int node, double newGroup, int layer, long numMCSteps, int numGroups);
std::vector<std::complex<double>> coeffs(int n);
MathMatrix bazziSample(size_t numNodes, size_t numLayers, size_t numGroups);
MathMatrix novelSample(size_t numNodes, size_t numLayers, size_t numGroups);
MathVector randomPartition(size_t n, int S);
long double lognewcorrectionNovel(MathMatrix G, int node, double newGroup, int layer, int numGroups, std::vector<std::vector<double>> corrections);
long double logNewProbNovel(MathMatrix G, int layer, int numGroups, std::vector<std::vector<double>> corrections);
long double logNewCorrection(MathVector prevGroups, MathVector nextGroups, int numGroups, std::vector<std::vector<double>> corrections);
long double logCorrectionBazziFull(MathMatrix G, MathVector newGroups, int layer, long numMCSteps, int numGroups);
long double logPDF(std::vector<MathMatrix> A, MathMatrix g, int layer);
double logfactorial2(double n);

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
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double perturbProb = 3e-3;
    std::random_device rd2;
    std::mt19937 gen2(rd2());
    std::uniform_int_distribution<> distrib2(0,numGroups-1);
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::uniform_int_distribution<> distrib4(0,numNodes-1);
    long numMCSteps = 1e3;
    std::vector<std::vector<double>> corrs = corrections(numNodes);
    
    std::vector<std::string> stepTypes = {"novel", "unif", "bazzi"};
    for (size_t sigmaIndex = 0; sigmaIndex < sigmaList.size(); ++ sigmaIndex)
    {
        double sigma = sigmaList[sigmaIndex];
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0,sigma);
        for (size_t group1Index = 0; group1Index < numRuns; ++group1Index)
        {
            std::vector<MathMatrix> A(layersA,MathMatrix(sizeA,sizeA));
            size_t group1Size = group1Sizes[group1Index];
            size_t group0Size = numNodes-group1Size;
            MathMatrix sbmProbs = {{sigmaList[sigmaIndex],0.10},{0.10,sigmaList[sigmaIndex]}};
            int groupI = 0;
            int groupJ = 0;
            
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
                std::ofstream myfile3("layer" + std::to_string(l) + "group1size" + std::to_string(group1Size) + "sigma" + std::to_string(sigmaList[sigmaIndex]) + "v2.txt",std::ios::trunc);
                myfile3 << A[l];
                myfile3.close();
            }
            for (size_t stepType = 0; stepType < 3; stepType++)
            {
                std::vector<MathVector> numCorrect(numLayers,MathVector(runs));
                std::vector<MathMatrix> initialGroups(numLayers,MathMatrix(numNodes,runs));
                std::vector<MathMatrix> finalGroups(numLayers,MathMatrix(numNodes,runs));
                
                for (size_t run = 0; run < runs; ++run)
                {
                    MathMatrix G = MathMatrix(numNodes,numLayers);
                    
                    if (stepTypes[stepType] == "unif")
                    {
                        for (size_t i = 0; i < numNodes; ++i)
                        {
                            for (size_t l = 0; l < numLayers; ++l)
                            {
                                G(i,l) = distrib2(gen2);
                            }
                        }
                    }
                    else if (stepTypes[stepType] == "bazzi")
                    {
                        G = bazziSample(numNodes, numLayers, numGroups);
                    }
                    else if (stepTypes[stepType] == "novel")
                    {
                        G = novelSample(numNodes, numLayers, numGroups);
                    }
                    for (size_t i = 0; i < numNodes; ++i)
                    {
                        for (size_t l = 0; l < numLayers; ++l)
                        {
                            initialGroups[l](i,run) = G(i,l);
                        }
                    }
                    for (size_t gibbsStep = 0; gibbsStep < 25; ++gibbsStep)
                    {
                        for (size_t i = 0; i < numNodes; ++i)
                        {
                            for (int l = 0; l < numLayers; ++l)
                            {
                                if (dis3(gen3) < perturbProb)
                                {
                                    if (l > 0)
                                    {
                                        double logproposalProb = 0.0;
                                        MathMatrix Gtemp = G;
                                        for (size_t lp = l; lp < numLayers; ++lp)
                                        {
                                            for (size_t k = 0; k < numNodes; ++k)
                                            {
                                                Gtemp(k,lp) = fabs(G(k,lp)-1.0);
                                            }
                                        }
                                        for (size_t lp = l; lp < numLayers; ++lp)
                                        {
                                            if (stepTypes[stepType] == "bazzi")
                                            {
                                                MathVector GtempL(numNodes);
                                                for (size_t k = 0; k < numNodes; ++k)
                                                {
                                                    GtempL(k) = Gtemp(k,lp);
                                                }
                                                logproposalProb += logCorrectionBazziFull(G, GtempL, lp, numMCSteps, numGroups);
                                            }
                                            else if (stepTypes[stepType] == "novel")
                                            {
                                                logproposalProb += logNewProbNovel(Gtemp, lp, numGroups, corrs) - logNewProbNovel(G, lp, numGroups, corrs);
                                            }
                                        }
                                        if (dis3(gen3) < std::exp(logproposalProb))
                                        {
                                            G = Gtemp;
                                        }
                                    }
                                }
                                else
                                {
                                    std::vector<double> logprobs(numGroups);
                                    std::vector<double> probs(numGroups);
                                    MathMatrix Gtemp = G;
                                    for (size_t g = 0; g < numGroups; ++g)
                                    {
                                        Gtemp(i,l) = g;
                                        logprobs[g] = logPDF(A, Gtemp, l);
                                        if (stepTypes[stepType] == "bazzi")
                                        {
                                            logprobs[g] += logCorrectionBazzi(G, i, g, l, numMCSteps, numGroups);
                                        }
                                        else if (stepTypes[stepType] == "novel")
                                        {
                                            logprobs[g] += lognewcorrectionNovel(G, i, g, l, numGroups, corrs);
                                        }
                                    }
                                    double maxlogProb = *max_element(logprobs.begin(), logprobs.end());
                                    double sum = 0;
                                    for (size_t g = 0; g < numGroups; ++g)
                                    {
                                        logprobs[g] -= maxlogProb;
                                        probs[g] = std::exp(logprobs[g]);
                                        sum += probs[g];
                                    }
                                    for (size_t g = 0; g < numGroups; ++g)
                                    {
                                        probs[g] /= sum;
                                    }
                                    double randomVal = dis(gen);
                                    sum = 0;
                                    int g = 0;
                                    bool stoppingCond = false;
                                    while (!stoppingCond)
                                    {
                                        if (randomVal < probs[g] + sum)
                                        {
                                            stoppingCond = true;
                                            G(i,l) = g;
                                        }
                                        else
                                        {
                                            sum += probs[g];
                                        }
                                        g++;
                                    }
                                }
                            }
                        }
                    }
                    std::cout << stepTypes[stepType] << " run " << run << " complete" << "\n";
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
                        std::ofstream myfile3("layer" + std::to_string(l) + "numCorrect" + std::to_string(group1Size) + stepTypes[stepType] + "sigma" + std::to_string(sigmaList[sigmaIndex]) +  "v2.txt",std::ios::trunc);
                        myfile3 << numCorrect[l];
                        myfile3.close();
                        std::ofstream myfile4("layer" + std::to_string(l) + "initialGroups" + std::to_string(group1Size) + stepTypes[stepType] + "sigma" + std::to_string(sigmaList[sigmaIndex]) + "v2.txt",std::ios::trunc);
                        myfile4 << initialGroups[l];
                        myfile4.close();
                        std::ofstream myfile5("layer" + std::to_string(l) + "finalGroups" + std::to_string(group1Size) + stepTypes[stepType] + "sigma" + std::to_string(sigmaList[sigmaIndex]) + "v2.txt",std::ios::trunc);
                        myfile5 << finalGroups[l];
                        myfile5.close();
                    }
                }
            }
        }
    }
    return 0;
}


long double logmcIntegrand(MathVector kronDel1, MathVector kronDel2, double p1, MathVector p2)
{
    size_t kdSize = kronDel1.size();
    long double out = 0.0;
    for (size_t i = 0; i < kdSize; ++i)
    {
        out += std::log(kronDel1(i)*p1 + (1.0-p1)*(p2(static_cast<int>(kronDel2(i)))));
    }
    return out;
}
long double sumlog(long double x, long double y)
{
    if (x > y)
    {
        return x + std::exp(y - x);
    }
    else
    {
        return y + std::exp(x - y);
    }
}
long double logbinom(double n, double k)
{
    return logfactorial2(n) - logfactorial2(k) - logfactorial2(n-k);
}
long double logCorrectionBazzi(MathMatrix G, int node, double newGroup, int layer, long numMCSteps, int numGroups)
{
    size_t numNodes = G.getRowSize();
    size_t numLayers = G.getColSize();
    long double logcorrectionnumprev = 0;
    long double logcorrectiondenomprev = 0;
    long double logcorrectionnumnext = 0;
    long double logcorrectiondenomnext = 0;
    MathVector krondel1numprev = MathVector(numNodes);
    MathVector krondel2numprev = MathVector(numNodes);
    MathVector krondel1denomprev = MathVector(numNodes);
    MathVector krondel2denomprev = MathVector(numNodes);
    MathVector krondel1numnext = MathVector(numNodes);
    MathVector krondel2numnext = MathVector(numNodes);
    MathVector krondel1denomnext = MathVector(numNodes);
    MathVector krondel2denomnext = MathVector(numNodes);
    for (size_t i = 0; i < numNodes; ++i)
    {
        if (layer > 0)
        {
            if (i == node)
            {
                krondel1numprev(i) = static_cast<double>(fabs(newGroup - G(i,layer-1)) < 0.1);
                krondel2numprev(i) = static_cast<double>(newGroup);
                krondel1denomprev(i) = static_cast<double>(fabs(G(i,layer) - G(i,layer-1)) < 0.1);
                krondel2denomprev(i) = static_cast<double>(G(i,layer));
            }
            else
            {
                krondel1numprev(i) = static_cast<double>(fabs(G(i,layer) - G(i,layer-1)) < 0.1);
                krondel2numprev(i) = static_cast<double>(G(i,layer));
                krondel1denomprev(i) = static_cast<double>(fabs(G(i,layer) - G(i,layer-1)) < 0.1);
                krondel2denomprev(i) = static_cast<double>(G(i,layer));
            }
        }
        if (layer < numLayers-1)
        {
            if (i == node)
            {
                krondel1numnext(i) = static_cast<double>(fabs(G(i,layer+1) - newGroup) < 0.1);
                krondel2numnext(i) = static_cast<double>(G(i,layer+1));
                krondel1denomnext(i) = static_cast<double>(fabs(G(i,layer+1) - G(i,layer)) < 0.1);
                krondel2denomnext(i) = static_cast<double>(G(i,layer+1));
            }
            else
            {
                krondel1numnext(i) = static_cast<double>(fabs(G(i,layer+1) - G(i,layer)) < 0.1);
                krondel2numnext(i) = static_cast<double>(G(i,layer+1));
                krondel1denomnext(i) = static_cast<double>(fabs(G(i,layer+1) - G(i,layer)) < 0.1);
                krondel2denomnext(i) = static_cast<double>(G(i,layer+1));
            }
        }
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::gamma_distribution<> dis4(1.0,1.0);
    for (size_t i = 0; i < numMCSteps; ++i)
    {
        double p1SampleNumprev = dis(gen);
        MathVector p2SampleNumprev = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleNumprev(j) = dis4(gen4);
        }
        p2SampleNumprev /= p2SampleNumprev.norm1();
        double p1SampleNumnext = dis(gen);
        MathVector p2SampleNumnext = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleNumnext(j) = dis4(gen4);
        }
        p2SampleNumnext /= p2SampleNumnext.norm1();
        double p1SampleDenomprev = dis(gen);
        MathVector p2SampleDenomprev = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleDenomprev(j) = dis4(gen4);
        }
        p2SampleDenomprev /= p2SampleDenomprev.norm1();
        double p1SampleDenomnext = dis(gen);
        MathVector p2SampleDenomnext = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleDenomnext(j) = dis4(gen4);
        }
        p2SampleDenomnext /= p2SampleDenomnext.norm1();
        if (i == 0)
        {
            if ((layer > 0) && (layer < numLayers-1))
            {
                logcorrectionnumprev = logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev);
                logcorrectiondenomprev = logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev);
                logcorrectionnumnext = logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext);
                logcorrectiondenomnext = logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext);
            }
            else if (layer <= 0)
            {
                logcorrectionnumnext = logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext);
                logcorrectiondenomnext = logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext);
            }
            else
            {
                logcorrectionnumprev = logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev);
                logcorrectiondenomprev = logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev);
            }
        }
        else
        {
            if ((layer > 0) && (layer < numLayers-1))
            {
                long double logcorrectionnumprevtemp = logcorrectionnumprev;
                logcorrectionnumprev = sumlog(logcorrectionnumprevtemp,logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev));
                long double logcorrectiondenomprevtemp = logcorrectiondenomprev;
                logcorrectiondenomprev = sumlog(logcorrectiondenomprevtemp,logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev));
                long double logcorrectionnumnexttemp = logcorrectionnumnext;
                logcorrectionnumnext = sumlog(logcorrectionnumnexttemp,logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext));
                long double logcorrectiondenomnexttemp = logcorrectiondenomnext;
                logcorrectiondenomnext = sumlog(logcorrectiondenomnexttemp, logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext));
            }
            else if (layer <= 0)
            {
                long double logcorrectionnumnexttemp = logcorrectionnumnext;
                logcorrectionnumnext = sumlog(logcorrectionnumnexttemp,logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext));
                long double logcorrectiondenomnexttemp = logcorrectiondenomnext;
                logcorrectiondenomnext = sumlog(logcorrectiondenomnexttemp, logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext));
            }
            else
            {
                long double logcorrectionnumprevtemp = logcorrectionnumprev;
                logcorrectionnumprev = sumlog(logcorrectionnumprevtemp,logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev));
                long double logcorrectiondenomprevtemp = logcorrectiondenomprev;
                logcorrectiondenomprev = sumlog(logcorrectiondenomprevtemp,logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev));
            }
        }
    }
    long double logcorrection = 0.0;
    if ((layer > 0) && (layer < numLayers-1))
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev) + (logcorrectionnumnext-logcorrectiondenomnext);
    }
    else if (layer <= 0)
    {
        logcorrection = (logcorrectionnumnext-logcorrectiondenomnext);
    }
    else
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev);
    }
    if (layer == 0)
    {
        double oldGroup = G(node,layer);
        double sizeoldGroup = 0.0;
        double sizenewGroup = 0.0;
        for (size_t j = 0; j < numNodes; ++j)
        {
            if (fabs(G(j,layer) - oldGroup) < 0.1)
            {
                sizeoldGroup += 1.0;
            }
            if (fabs(G(j,layer) - newGroup) < 0.1)
            {
                sizenewGroup += 1.0;
            }
        }
        logcorrection += std::log(sizenewGroup + 1) - std::log(sizeoldGroup);
    }
    return logcorrection;
}

long double logCorrectionBazziFull(MathMatrix G, MathVector newGroups, int layer, long numMCSteps, int numGroups)
{
    size_t numNodes = G.getRowSize();
    size_t numLayers = G.getColSize();
    long double logcorrectionnumprev = 0;
    long double logcorrectiondenomprev = 0;
    long double logcorrectionnumnext = 0;
    long double logcorrectiondenomnext = 0;
    MathVector krondel1numprev = MathVector(numNodes);
    MathVector krondel2numprev = MathVector(numNodes);
    MathVector krondel1denomprev = MathVector(numNodes);
    MathVector krondel2denomprev = MathVector(numNodes);
    MathVector krondel1numnext = MathVector(numNodes);
    MathVector krondel2numnext = MathVector(numNodes);
    MathVector krondel1denomnext = MathVector(numNodes);
    MathVector krondel2denomnext = MathVector(numNodes);
    for (size_t i = 0; i < numNodes; ++i)
    {
        if (layer > 0)
        {
            krondel1numprev(i) = static_cast<double>(fabs(newGroups(i) - G(i,layer-1)) < 0.1);
            krondel2numprev(i) = static_cast<double>(newGroups(i));
            krondel1denomprev(i) = static_cast<double>(fabs(G(i,layer) - G(i,layer-1)) < 0.1);
            krondel2denomprev(i) = static_cast<double>(G(i,layer));
            
        }
        if (layer < numLayers-1)
        {
            krondel1numnext(i) = static_cast<double>(fabs(G(i,layer+1) - newGroups(i)) < 0.1);
            krondel2numnext(i) = static_cast<double>(G(i,layer+1));
            krondel1denomnext(i) = static_cast<double>(fabs(G(i,layer+1) - G(i,layer)) < 0.1);
            krondel2denomnext(i) = static_cast<double>(G(i,layer+1));
        }
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::gamma_distribution<> dis4(1.0,1.0);
    for (size_t i = 0; i < numMCSteps; ++i)
    {
        double p1SampleNumprev = dis(gen);
        MathVector p2SampleNumprev = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleNumprev(j) = dis4(gen4);
        }
        p2SampleNumprev /= p2SampleNumprev.norm1();
        double p1SampleNumnext = dis(gen);
        MathVector p2SampleNumnext = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleNumnext(j) = dis4(gen4);
        }
        p2SampleNumnext /= p2SampleNumnext.norm1();
        double p1SampleDenomprev = dis(gen);
        MathVector p2SampleDenomprev = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleDenomprev(j) = dis4(gen4);
        }
        p2SampleDenomprev /= p2SampleDenomprev.norm1();
        double p1SampleDenomnext = dis(gen);
        MathVector p2SampleDenomnext = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p2SampleDenomnext(j) = dis4(gen4);
        }
        p2SampleDenomnext /= p2SampleDenomnext.norm1();
        if (i == 0)
        {
            if ((layer > 0) && (layer < numLayers-1))
            {
                logcorrectionnumprev = logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev);
                logcorrectiondenomprev = logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev);
                logcorrectionnumnext = logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext);
                logcorrectiondenomnext = logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext);
            }
            else if (layer <= 0)
            {
                logcorrectionnumnext = logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext);
                logcorrectiondenomnext = logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext);
            }
            else
            {
                logcorrectionnumprev = logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev);
                logcorrectiondenomprev = logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev);
            }
        }
        else
        {
            if ((layer > 0) && (layer < numLayers-1))
            {
                long double logcorrectionnumprevtemp = logcorrectionnumprev;
                logcorrectionnumprev = sumlog(logcorrectionnumprevtemp,logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev));
                long double logcorrectiondenomprevtemp = logcorrectiondenomprev;
                logcorrectiondenomprev = sumlog(logcorrectiondenomprevtemp,logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev));
                long double logcorrectionnumnexttemp = logcorrectionnumnext;
                logcorrectionnumnext = sumlog(logcorrectionnumnexttemp,logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext));
                long double logcorrectiondenomnexttemp = logcorrectiondenomnext;
                logcorrectiondenomnext = sumlog(logcorrectiondenomnexttemp, logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext));
            }
            else if (layer <= 0)
            {
                long double logcorrectionnumnexttemp = logcorrectionnumnext;
                logcorrectionnumnext = sumlog(logcorrectionnumnexttemp,logmcIntegrand(krondel1numnext, krondel2numnext, p1SampleNumnext, p2SampleNumnext));
                long double logcorrectiondenomnexttemp = logcorrectiondenomnext;
                logcorrectiondenomnext = sumlog(logcorrectiondenomnexttemp, logmcIntegrand(krondel1denomnext, krondel2denomnext, p1SampleDenomnext, p2SampleDenomnext));
            }
            else
            {
                long double logcorrectionnumprevtemp = logcorrectionnumprev;
                logcorrectionnumprev = sumlog(logcorrectionnumprevtemp,logmcIntegrand(krondel1numprev, krondel2numprev, p1SampleNumprev, p2SampleNumprev));
                long double logcorrectiondenomprevtemp = logcorrectiondenomprev;
                logcorrectiondenomprev = sumlog(logcorrectiondenomprevtemp,logmcIntegrand(krondel1denomprev, krondel2denomprev, p1SampleDenomprev, p2SampleDenomprev));
            }
        }
    }
    long double logcorrection = 0.0;
    if ((layer > 0) && (layer < numLayers-1))
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev);
    }
    else if (layer <= 0)
    {
        logcorrection = 0;
    }
    else
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev);
    }
    return logcorrection;
}

long double logPDF(std::vector<MathMatrix> A, MathMatrix g, int layer)
{
    long double out = 0.0;
    MathMatrix M(2, 2, 0.0);
    MathMatrix T(2, 2, 0.0);
    size_t numNodes = A[layer].getColSize();
    for (size_t i = 0; i < numNodes; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            T(g(i,layer),g(j,layer)) += 1.0;
            M(g(i,layer),g(j,layer)) += A[layer](i,j);
        }
    }
    out += logfactorial2(M(0,0)) + logfactorial2(T(0,0) - M(0,0)) - logfactorial2(T(0,0) + 1);
    out += logfactorial2(M(0,1) + M(1,0)) + logfactorial2(T(0,1) + T(1,0) - M(0,1) - M(1,0)) - logfactorial2(T(0,1) + T(1,0) + 1);
    out += logfactorial2(M(1,1)) + logfactorial2(T(1,1) - M(1,1)) - logfactorial2(T(1,1) + 1);
    return out;
}

std::vector<std::complex<double>> coeffs(int n)
{
    std::vector<std::complex<double>> out(n+1,{0,0});
    std::complex<double> last = {1,0};
    for (size_t k = 0; k <= n-1; ++k)
    {
        std::complex<double> sum = {0,0};
        for (size_t r = 1; r <= n; ++r)
        {
            std::complex<double> temp = {1,0};
            for (size_t m = 1; m <= n; ++m)
            {
                if (m != r)
                {
                    temp /= std::polar(1.0, 2.0 * M_PI * r / (n+1)) - std::polar(1.0, 2.0 * M_PI * m / (n+1));
                }
            }
            sum += 1.0 * temp * std::polar(1.0,2.0 * M_PI * r * k / (n+1)) * (std::log(1.0 - std::polar(1.0,2.0 * M_PI * r  / (n+1))) - std::log(-1.0 * std::polar(1.0, 2.0 * M_PI * r  / (n+1))));
        }
        out[k] = sum;
        last -= out[k];
    }
    out[n] = last;
    return out;
}

std::vector<std::vector<double>> corrections(int n)
{
    std::vector<double> outTemp(n+1,0);
    std::vector<std::vector<double>> out(n+1,outTemp);
    for (size_t k = 1; k <= n; ++k)
    {
        bool stoppingCond = false;
        std::vector<std::complex<double>> coeffsN = coeffs(k);
        for (size_t a = 0; a < k + 1; ++a)
        {
            if (!stoppingCond)
            {
                out[k][a] = std::log(coeffsN[a].real());
            }
            else
            {
                out[k][a] = out[k][a-1] - std::log(1.0);
            }
            stoppingCond = (out[k][a] < -16.0);
        }
    }
    return out;
}

long double logNewCorrection(MathVector prevGroups, MathVector nextGroups, int numGroups, std::vector<std::vector<double>> corrections)
{
    long double out = 0.0;
    size_t numNodes = prevGroups.size();
    for (size_t i = 0; i < numGroups; ++i)
    {
        int groupSize = 0;
        int numSame = 0;
        MathVector groupSizes = MathVector(numGroups-1);
        for (size_t j = 0; j < numNodes; ++j)
        {
            if (fabs(prevGroups(j) - i) < 0.1)
            {
                groupSize++;
                if (fabs(nextGroups(j) - i) < 0.1)
                {
                    numSame++;
                }
                else if (nextGroups(j) < i)
                {
                    groupSizes((int)nextGroups(j)) = groupSizes((int)nextGroups(j)) + 1.0;
                }
                else
                {
                    groupSizes((int)nextGroups(j)-1) = groupSizes((int)nextGroups(j)-1) + 1.0;
                }
            }
        }
        out += corrections[groupSize][groupSize-numSame] - logbinom(groupSize-numSame+numGroups-2, numGroups-2) - logfactorial2(groupSize);
        for (size_t k = 0; k < numGroups-1; ++k)
        {
            out += logfactorial2(groupSizes(k));
        }
        out += logfactorial2(numSame);
    }
    return out;
}
MathMatrix bazziSample(size_t numNodes, size_t numLayers, size_t numGroups)
{
    std::random_device rd3;
    std::mt19937 gen3(rd3());
    std::uniform_real_distribution<> dis3(0.0, 1.0);
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::gamma_distribution<> dis4(1.0,1.0);
    MathMatrix G = MathMatrix(numNodes,numLayers);
    std::vector<double> firstLayerTemp(numNodes);
    MathVector firstLayerPartition = randomPartition(numGroups, numNodes);
    int offset = 0;
    for (size_t k = 0; k < numGroups; ++k)
    {
        for (size_t i = 0; i < firstLayerPartition(k); ++i)
        {
            firstLayerTemp[i + offset] = k;
        }
        offset += firstLayerPartition(k);
    }
    std::vector<int> randPerm;
    for (int i=0; i<numNodes; ++i)
    {
        randPerm.push_back(i);
    }
    std::random_shuffle(randPerm.begin(), randPerm.end());
    for (int i=0; i<numNodes; ++i)
    {
        G(i,0) = firstLayerTemp[randPerm[i]];
    }
    
    for (size_t l = 1; l < numLayers; ++l)
    {
        double P = dis3(gen3);
        MathVector p = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p(j) = dis4(gen4);
        }
        p /= p.norm1();
        for (size_t i = 0; i < numNodes; ++i)
        {
            if (dis3(gen3) < P)
            {
                G(i,l) = G(i,l-1);
            }
            else
            {
                double sample = dis3(gen3);
                double group = 0;
                double tempProb = 0;
                for (size_t m = 0; m < numGroups; ++m)
                {
                    tempProb += p(m);
                    if (sample > tempProb)
                    {
                        group += 1;
                    }
                }
                G(i,l) = group;
            }
        }
    }
    return G;
}
MathMatrix novelSample(size_t numNodes, size_t numLayers, size_t numGroups)
{
    std::random_device rd3;
    std::mt19937 gen3(rd3());
    std::uniform_real_distribution<> dis3(0.0, 1.0);
    std::random_device rd4;
    std::mt19937 gen4(rd4());
    std::gamma_distribution<> dis4(1.0,1.0);
    MathMatrix G = MathMatrix(numNodes,numLayers);
    std::vector<double> firstLayerTemp(numNodes);
    MathVector firstLayerPartition = randomPartition(numGroups, numNodes);
    int offset = 0;
    for (size_t k = 0; k < numGroups; ++k)
    {
        for (size_t i = 0; i < firstLayerPartition(k); ++i)
        {
            firstLayerTemp[i + offset] = k;
        }
        offset += firstLayerPartition(k);
    }
    std::vector<int> randPerm;
    for (int i=0; i<numNodes; ++i)
    {
        randPerm.push_back(i);
    }
    std::random_shuffle(randPerm.begin(), randPerm.end());
    for (int i=0; i<numNodes; ++i)
    {
        G(i,0) = firstLayerTemp[randPerm[i]];
    }
    for (size_t l = 1; l < numLayers; ++l)
    {
        for (size_t k = 0; k < numGroups; ++k)
        {
            double P = dis3(gen3);
            std::vector<int> indices;
            int groupSize = 0;
            for (size_t i = 0; i < numNodes; ++i)
            {
                if (fabs(G(i,l-1)-k) < 0.1)
                {
                    indices.push_back(i);
                    groupSize++;
                }
            }
            MathVector p = MathVector(groupSize+1);
            for (size_t j = 0; j < groupSize+1; ++j)
            {
                p(j) = std::pow(P,groupSize-j);
            }
            p /= p.norm1();
            double sample = dis3(gen3);
            int numSame = 0;
            double tempProb = 0;
            for (size_t m = 0; m < groupSize+1; ++m)
            {
                tempProb += p(m);
                if (sample > tempProb)
                {
                    numSame += 1;
                }
            }
            MathVector groupPartition = randomPartition(numGroups-1, groupSize-numSame);
            std::vector<double> groupTemp;
            for (size_t g = 0; g < numGroups; ++g)
            {
                if (g < k)
                {
                    for (size_t i = 0; i < groupPartition(g); ++i)
                    {
                        groupTemp.push_back(g);
                    }
                }
                else if (g > k)
                {
                    for (size_t i = 0; i < groupPartition(g-1); ++i)
                    {
                        groupTemp.push_back(g);
                    }
                }
                else
                {
                    for (size_t i = 0; i < numSame; ++i)
                    {
                        groupTemp.push_back(k);
                    }
                }
            }
            std::vector<int> randGroupPerm;
            for (int i=0; i<groupSize; ++i)
            {
                randGroupPerm.push_back(i);
            }
            std::random_shuffle(randGroupPerm.begin(), randGroupPerm.end());
            for (int i=0; i<groupSize; ++i)
            {
                G(indices[i],l) = groupTemp[randGroupPerm[i]];
            }
        }
        double P = dis3(gen3);
        MathVector p = MathVector(numGroups);
        for (size_t j = 0; j < numGroups; ++j)
        {
            p(j) = dis4(gen4);
        }
        p /= p.norm1();
        for (size_t i = 0; i < numNodes; ++i)
        {
            if (dis3(gen3) < P)
            {
                G(i,l) = G(i,l-1);
            }
            else
            {
                double sample = dis3(gen3);
                double group = 0;
                double tempProb = 0;
                for (size_t m = 0; m < numGroups; ++m)
                {
                    tempProb += p(m);
                    if (sample > tempProb)
                    {
                        group += 1;
                    }
                }
                G(i,l) = group;
            }
        }
    }
    return G;
}
MathVector randomPartition(size_t n, int S)
{
    MathMatrix P = MathMatrix(S+1,n,1.0);
    std::random_device rd3;
    std::mt19937 gen3(rd3());
    std::uniform_real_distribution<> dis3(0.0, 1.0);
    for (size_t i = 0; i < S+1; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            for (size_t k = 1; k < i+1; ++k)
            {
                P(i,j) = P(i,j) * (n+i-j-k) / k;
            }
        }
    }
    int s = S;
    MathVector R = MathVector(n);
    for (size_t j = 0; j < n; ++j)
    {
        int count = 0;
        double threshold = dis3(gen3) * P(s,j);
        for (size_t k = 0; k < s; ++k)
        {
            count += (int)(threshold <= P(k,j));
        }
        R(j) = count;
        s -= count;
    }
    return R;
}
long double lognewcorrectionNovel(MathMatrix G, int node, double newGroup, int layer, int numGroups, std::vector<std::vector<double>> corrections)
{
    size_t numNodes = G.getRowSize();
    size_t numLayers = G.getColSize();
    long double logcorrectionnumprev = 0;
    long double logcorrectiondenomprev = 0;
    long double logcorrectionnumnext = 0;
    long double logcorrectiondenomnext = 0;
    MathVector prevLayer = MathVector(numNodes);
    MathVector currLayerOld = MathVector(numNodes);
    MathVector currLayerNew = MathVector(numNodes);
    MathVector nextLayer = MathVector(numNodes);
    if (layer > 0)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            prevLayer(i) = G(i,layer-1);
        }
    }
    if (layer < numLayers-1)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            nextLayer(i) = G(i,layer+1);
        }
    }
    for (size_t i = 0; i < numNodes; ++i)
    {
        if (i == node)
        {
            currLayerOld(i) = G(i,layer);
            currLayerNew(i) = newGroup;
        }
        else
        {
            currLayerOld(i) = G(i,layer);
            currLayerNew(i) = G(i,layer);
        }
    }
    if (layer > 0)
    {
        logcorrectionnumprev = logNewCorrection(prevLayer, currLayerNew, numGroups, corrections);
        logcorrectiondenomprev = logNewCorrection(prevLayer, currLayerOld, numGroups, corrections);
    }
    if (layer < numLayers-1)
    {
        logcorrectionnumnext = logNewCorrection(currLayerNew, nextLayer, numGroups, corrections);
        logcorrectiondenomnext = logNewCorrection(currLayerOld, nextLayer, numGroups, corrections);
    }
    long double logcorrection = 0.0;
    if ((layer > 0) && (layer < numLayers-1))
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev) + (logcorrectionnumnext-logcorrectiondenomnext);
    }
    else if (layer <= 0)
    {
        logcorrection = (logcorrectionnumnext-logcorrectiondenomnext);
    }
    else
    {
        logcorrection = (logcorrectionnumprev-logcorrectiondenomprev);
    }
    if (layer == 0)
    {
        double oldGroup = G(node,layer);
        double sizeoldGroup = 0.0;
        double sizenewGroup = 0.0;
        for (size_t j = 0; j < numNodes; ++j)
        {
            if (fabs(G(j,layer) - oldGroup) < 0.1)
            {
                sizeoldGroup += 1.0;
            }
            if (fabs(G(j,layer) - newGroup) < 0.1)
            {
                sizenewGroup += 1.0;
            }
        }
        logcorrection += std::log(sizenewGroup + 1) - std::log(sizeoldGroup);
    }
    return logcorrection;
}

long double logNewProbNovel(MathMatrix G, int layer, int numGroups, std::vector<std::vector<double>> corrections)
{
    size_t numNodes = G.getRowSize();
    size_t numLayers = G.getColSize();
    long double logcorrectiondenomprev = 0;
    long double logcorrectiondenomnext = 0;
    MathVector prevLayer = MathVector(numNodes);
    MathVector currLayerOld = MathVector(numNodes);
    MathVector nextLayer = MathVector(numNodes);
    if (layer > 0)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            prevLayer(i) = G(i,layer-1);
        }
    }
    if (layer < numLayers-1)
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            nextLayer(i) = G(i,layer+1);
        }
    }
    for (size_t i = 0; i < numNodes; ++i)
    {
        currLayerOld(i) = G(i,layer);
    }
    if (layer > 0)
    {
        logcorrectiondenomprev = logNewCorrection(prevLayer, currLayerOld, numGroups, corrections);
    }
    if (layer < numLayers-1)
    {
        logcorrectiondenomnext = logNewCorrection(currLayerOld, nextLayer, numGroups, corrections);
    }
    long double logcorrection = 0.0;
    if ((layer > 0) && (layer < numLayers-1))
    {
        logcorrection = logcorrectiondenomprev;
    }
    else if (layer <= 0)
    {
        logcorrection = 0;
    }
    else
    {
        logcorrection = logcorrectiondenomprev;
    }
    if (layer == 0)
    {
        MathVector groupSizes = MathVector(numGroups);
        for (size_t j = 0; j < numNodes; ++j)
        {
            groupSizes(static_cast<int>(G(j,layer))) = groupSizes(static_cast<int>(G(j,layer))) + 1.0;
        }
        for (size_t g = 0; g < numGroups; ++g)
        {
            logcorrection += logfactorial2(groupSizes(g));
        }
    }
    return logcorrection;
}

double logfactorial2(double n)
{
    long double out = 0.0;
    if (n > 1)
    {
        out = lgamma(1.0 + (long double)n);
    }
    return out;
}

