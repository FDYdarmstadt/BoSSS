GenerateCubeFiles
generates cube-files for further usage with Extra-P
ImbalanceDetection
Computes median and average of cpuwise runtimes of a method based on method call record (using ilPSP.Tracing). This can be used as an first indicator of load imbalances over cores. Note that, if median and average are not similar, this is a sufficiant criterion for imbalance. But average similar median may not ensure, that there is no imblance.
GenerateRuntimeCluster
clustering of runtime of specific methods via method call record (using ilPSP.Tracing). 2 available algorithms (k-Means and relative gap). The relative gap algorithm clusters data into N cluster of an 1D space according to the Nth biggest jumps of data values. Means: data values are sorted in ascending order. (data(x+1)-data(x))/